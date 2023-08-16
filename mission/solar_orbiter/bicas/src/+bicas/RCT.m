%
% Class that collects generic functions related to reading different kinds of
% RCTs (files) into corresponding data structures without any real modification.
%
%
% DESIGN INTENT
% =============
% Implemented so that no calibration data is modified/added to/removed from. The
% returned data structures reflect the content of the RCTs, not necessarily the
% data used. Modification of data should be done elsewhere, in particular
% modifications of transfer functions, e.g. extrapolation, cut-offs,
% inversions.
% --
% NOTE: BIAS & LFR RCTs: contain FTFs which are not inverted in this code.
%       TDS RCTs:        contain ITFs.
% NOTE: Code still converts RCT TFs slightly:
%   frequency      : Hz    --> rad/s
%   phase+amplitude: degrees,dimensionless real value --> Z (complex number)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-15
%
classdef RCT
% PROPOSAL: Rename.
%   CON: Does not conform to naming conventions.
%   PRO: Only contains functions for reading RCTs (but not modifying).
%        Cf bicas.proc.L1L2.cal_RCT_types.
%       CON: Irrelevant. All RCT-code is about reading RCTs, kind of.
%   PROPOSAL: rctread
%
% PROPOSAL: Merge with bicas.proc.L1L2.cal_RCT_types.
% PROPOSAL: Class for each type of RCT. Common superclass. Can contain 
%   (1) methods read_*_RCT()
%   (2) methods modify_*_data()
%   (3) methods log_*_RCTs()
%   (4) filenameRegexpSettingKey
%   PRO: Can replace
%       (1) Structs created by
%           bicas.proc.L1L2.cal_RCT_types.init_RCT_TYPES_MAP.entry().
%       (2) Bulk of bicas.RCT, bicas.proc.L1L2.cal_RCT_types
%       (3) RctData structs returned by read_*_RCT() and modify_*_data().
%       CON: Might use shared private functions that need to live in some other
%            file.
%       CON: There is a need to be able to load the "plain" content of RCTs,
%            without modification.
%           Ex: Best practice w.r.t. modularization. Reusability.
%           Ex: Loading and analyzing RCTs outside of BICAS.
%           PROPOSAL: Can still have separate methods for (1) reading raw RCT, and (2)
%                     modifying the raw RCT data.
%               CON: Still one class for two datastructures (per RCT type).
%                   CON-PROPOSAL: Store both raw AND modified RCT data in class.
%                   CON-PROPOSAL: Store only raw OR  modified RCT data in class.
%
% PROPOSAL: Move bicas.RCT
%   CON: Contains generic RCT functionality. Not directly processing related.
%     --> bicas.proc.L1L2.RCT ?
%     --> bicas.proc.L1L2.cal.RCT ?
%     --> bicas.proc.L1L2*.RCT_read ?
%     --> bicas.RCT_read ?
%     --> bicas.read_RCT ?
%   PRO: RCTs are specific to L1/L1R-->L2 processing.
%    
% PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Assert CDF skeleton/master version number.
% PROPOSAL: Assert skeleton/master.
%   PRO: Can give better error message when reading the wrong RCT.
%
% PROPOSAL: Function for permuting indices to handle dataobj's handling of 1 record-case.
%
% PROPOSAL: Assert/warn (depending on setting?) when CDF metadata imply that the RCT zVariables have the wrong units.
% PROPOSAL: Use utility function for reading every zVariable.
%   PROPOSAL: Assert units from zVar attributes.
%
% PROPOSAL: Classes for RCT data.
%   PRO: BIAS data has many fields.
%   PRO: More well-defined data structs.
%   PRO: Automatic assertions.
%   CON: Structs are modified when cal.m uses them, i.e. one could just as well
%        have classes for the format cal.m uses. ==> Too many classes.



    properties(Access=private, Constant)

        % Minimum number of numerator or denominator coefficients in the BIAS RCT.
        N_MIN_TF_NUMER_DENOM_COEFFS = 8;
        
        % Minimum number of expected entries in tabulated transfer functions in RCTs.
        TF_TABLE_MIN_LENGTH = 10;
        
    end


    
    methods(Static, Access=public)
        
        
        
        function [RctData] = read_BIAS_RCT(filePath)
            % TODO-DEC: How handle time?
            %   PROPOSAL: "Only" access the BIAS values (trans.func and other) through a function instead of selecting
            %             indices in a data struct.
            %       PROPOSAL: (private method) [omegaRps, zVpc] = get_transfer_func(epoch, signalType)
            %           signalType = 'DC single' etc
            
            Do = dataobj(filePath);
            
            % Constants for interpreting the array indices in the CDF.
            I_NUMERATOR   = 1;
            I_DENOMINATOR = 2;
            %
            I_DC_SINGLE = 1;
            I_DC_DIFF   = 2;
            I_AC_LG     = 3;
            I_AC_HG     = 4;
            %
            I_E12 = 1;
            I_E13 = 2;
            I_E23 = 3;
            
            try
                % NOTE: Assumes 1 CDF record or many (time-dependent values).
                % ==> Must handle that dataobj assigns differently for these two
                %     cases.
                epochL                    = bicas.RCT.normalize_dataobj_ZV(Do.data.Epoch_L);
                epochH                    = bicas.RCT.normalize_dataobj_ZV(Do.data.Epoch_H);
                biasCurrentOffsetsAAmpere = bicas.RCT.normalize_dataobj_ZV(Do.data.BIAS_CURRENT_OFFSET);      % DEPEND_0 = Epoch_L
                biasCurrentGainsAapt      = bicas.RCT.normalize_dataobj_ZV(Do.data.BIAS_CURRENT_GAIN);        % DEPEND_0 = Epoch_L
                dcSingleOffsetsAVolt      = bicas.RCT.normalize_dataobj_ZV(Do.data.V_OFFSET);                 % DEPEND_0 = Epoch_H
                dcDiffOffsetsAVolt        = bicas.RCT.normalize_dataobj_ZV(Do.data.E_OFFSET);                 % DEPEND_0 = Epoch_H
                ftfCoeffs                 = bicas.RCT.normalize_dataobj_ZV(Do.data.TRANSFER_FUNCTION_COEFFS); % DEPEND_0 = Epoch_L

                nEpochL = size(epochL, 1);
                nEpochH = size(epochH, 1);

                % IMPLEMENTATION NOTE: Corrects for what seems to be a bug in
                % dataobj. dataobj permutes/removes indices, and permutes them
                % differently depending on the number of CDF records (but wrong
                % in all cases).
                %
                % 1 CDF record : cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       F/TTT"   # 3=number of dimensions/record
                % 2 CDF records: cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       T/TTT"
                % 1 CDF record:  size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [4 2 8]
                % 2 CDF records: size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [2 4 2 8]
                ftfCoeffs = permute(ftfCoeffs, [1, 4, 3, 2]);



                %=======================================================
                % ASSERTIONS: Size of tfCoeffs/TRANSFER_FUNCTION_COEFFS
                %=======================================================
                % ND = Numerator Denominator
                nNdCoeffs = irf.assert.sizes(ftfCoeffs, [nEpochL, -1, 2, 4]);
                assert(nNdCoeffs >= bicas.RCT.N_MIN_TF_NUMER_DENOM_COEFFS)

                %================================
                % Assign struct that is returned
                %================================
                RctData.epochL = epochL;
                RctData.epochH = epochH;

                RctData.Current.offsetsAAmpere = biasCurrentOffsetsAAmpere;
                RctData.Current.gainsAapt      = biasCurrentGainsAapt;
                RctData.dcSingleOffsetsAVolt   = dcSingleOffsetsAVolt;
                RctData.DcDiffOffsets.E12AVolt = dcDiffOffsetsAVolt(:, I_E12);
                RctData.DcDiffOffsets.E13AVolt = dcDiffOffsetsAVolt(:, I_E13);
                RctData.DcDiffOffsets.E23AVolt = dcDiffOffsetsAVolt(:, I_E23);

                % NOTE: Using name "FtfSet" only to avoid "Ftfs" (plural).
                % (List, Table would be wrong? Use "FtfTable"?)
                RctData.FtfSet.DcSingleAvpiv = bicas.RCT.create_TF_sequence(...
                    ftfCoeffs(:, :, I_NUMERATOR,   I_DC_SINGLE), ...
                    ftfCoeffs(:, :, I_DENOMINATOR, I_DC_SINGLE));

                RctData.FtfSet.DcDiffAvpiv = bicas.RCT.create_TF_sequence(...
                    ftfCoeffs(:, :, I_NUMERATOR,   I_DC_DIFF), ...
                    ftfCoeffs(:, :, I_DENOMINATOR, I_DC_DIFF));

                RctData.FtfSet.AcLowGainAvpiv = bicas.RCT.create_TF_sequence(...
                    ftfCoeffs(:, :, I_NUMERATOR,   I_AC_LG), ...
                    ftfCoeffs(:, :, I_DENOMINATOR, I_AC_LG));

                RctData.FtfSet.AcHighGainAvpiv = bicas.RCT.create_TF_sequence(...
                    ftfCoeffs(:, :, I_NUMERATOR,   I_AC_HG), ...
                    ftfCoeffs(:, :, I_DENOMINATOR, I_AC_HG));
                
                % ASSERTIONS
                irf.assert.sizes(...
                    RctData.FtfSet.DcSingleAvpiv,   [nEpochL, 1], ...
                    RctData.FtfSet.DcDiffAvpiv,     [nEpochL, 1], ...
                    RctData.FtfSet.AcLowGainAvpiv,  [nEpochL, 1], ...
                    RctData.FtfSet.AcHighGainAvpiv, [nEpochL, 1]);
                for iEpochL = 1:nEpochL
                    %assert(Bias.ItfSet.DcSingleAvpiv{iEpochL}.eval(0) > 0, 'BICAS:FailedToReadInterpretRCT', 'DC single inverted transfer function is not positive (and real) at 0 Hz. (Wrong sign?)');
                    %assert(Bias.ItfSet.DcDiffAvpiv{iEpochL}.eval(0)   > 0, 'BICAS:FailedToReadInterpretRCT',   'DC diff inverted transfer function is not positive (and real) at 0 Hz. (Wrong sign?)');
                    % Unsure if assertion makes sense for AC, or possibly even
                    % for DC.
                    % 2020-03-10: This criterion is not true for AC high-gain
                    % transfer function fit now used (but does for AC diff
                    % low-gain).
                end
                
                %==============================================================
                % ASSERTIONS:
                % All variables NOT based on tfCoeffs/TRANSFER_FUNCTION_COEFFS
                %==============================================================
                bicas.utils.assert_ZV_Epoch(RctData.epochL)
                bicas.utils.assert_ZV_Epoch(RctData.epochH)
                validateattributes(RctData.epochL, {'numeric'}, {'increasing'})
                validateattributes(RctData.epochH, {'numeric'}, {'increasing'})

                irf.assert.sizes(...
                    RctData.Current.offsetsAAmpere, [nEpochL, 3], ...
                    RctData.Current.gainsAapt,      [nEpochL, 3], ...
                    RctData.dcSingleOffsetsAVolt,   [nEpochH, 3]);

                for fn = fieldnames(RctData.DcDiffOffsets)'
                    irf.assert.sizes(RctData.DcDiffOffsets.(fn{1}), [nEpochH, 1]);
                end
                
            catch Exc1
                Exc2 = MException(...
                    'BICAS:FailedToReadInterpretRCT', ...
                    'Error when interpreting calibration file (BIAS RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2)
            end
        end



        % LfrFtfTpivTable : {iLsf}{iBlts}. Table of LFR FTFs.
        %                   iLsf=1..3 : iBlts=1..5 for BLTS 1-5.
        %                   iLsf=4    : iBlts=1..3 for BIAS 1-3.
        function RctData = read_LFR_RCT(filePath)
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special
                % behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.

                % NOTE: There are separate TFs for each BLTS channel, not just
                % separate LFR sampling frequencies, i.e. there are 5+5+5+3=18
                % TFs (but only 1 frequency table/LSF, since they are recycled).
                % NOTE: The assignment of indices here effectively determines
                % the translation between array index and LFR Sampling Frequency
                % (LSF). This is NOT the same as the values in the LFR zVar
                % FREQ.
                freqTableHz{1}   = shiftdim(Do.data.Freqs_F0.data);    % NOTE: Index {iLsf}.
                freqTableHz{2}   = shiftdim(Do.data.Freqs_F1.data);
                freqTableHz{3}   = shiftdim(Do.data.Freqs_F2.data);
                freqTableHz{4}   = shiftdim(Do.data.Freqs_F3.data);

                amplTableTpiv{1} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F0.data);
                amplTableTpiv{2} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F1.data);
                amplTableTpiv{3} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F2.data);
                amplTableTpiv{4} = shiftdim(Do.data.TF_BIAS_123_amplitude_F3.data);

                phaseTableDeg{1} = shiftdim(Do.data.TF_BIAS_12345_phase_F0.data);
                phaseTableDeg{2} = shiftdim(Do.data.TF_BIAS_12345_phase_F1.data);
                phaseTableDeg{3} = shiftdim(Do.data.TF_BIAS_12345_phase_F2.data);
                phaseTableDeg{4} = shiftdim(Do.data.TF_BIAS_123_phase_F3.data);

                for iLsf = 1:4
                    % NOTE: F3 is an exception and has no AC (iBlts={4,5}) TF.
                    if iLsf ~= 4   nBlts = 5;
                    else           nBlts = 3;
                    end

                    % NOTE: Values for the specific LSF, hence the prefix.
                    lsfFreqTableHz   = freqTableHz{  iLsf};
                    lsfAmplTableTpiv = amplTableTpiv{iLsf};
                    lsfPhaseTableDeg = phaseTableDeg{iLsf};

                    % ASSERTIONS: Check CDF array sizes, and implicitly that the
                    % CDF format is the expected one.
                    nFreqs = irf.assert.sizes(...
                        lsfFreqTableHz,   [-1,       1 ], ...
                        lsfAmplTableTpiv, [-1, nBlts], ...
                        lsfPhaseTableDeg, [-1, nBlts]);
                    assert(nFreqs >= bicas.RCT.TF_TABLE_MIN_LENGTH)

                    for iBlts = 1:nBlts
                        
                        lsfBltsFreqTableHz   = lsfFreqTableHz;
                        lsfBltsAmplTableTpiv = lsfAmplTableTpiv(:, iBlts);
                        lsfBltsPhaseTableDeg = lsfPhaseTableDeg(:, iBlts);

                        FtfTpiv = irf.utils.tabulated_transform(...
                            lsfBltsFreqTableHz * 2*pi, ...
                            lsfBltsAmplTableTpiv, ...
                            deg2rad(lsfBltsPhaseTableDeg));
                        
                        % ASSERTION: FTF
                        assert(FtfTpiv.toward_zero_at_high_freq())
                        
                        FtfTpivTable{iLsf}{iBlts} = FtfTpiv;
                    end
                end
                
                % NOTE: Storing data in struct field to clarify the nature of
                % the content to the caller.
                RctData = [];
                RctData.FtfTpivTable = FtfTpivTable;
                
            catch Exc1
                Exc2 = MException(...
                    'BICAS:FailedToReadInterpretRCT', ...
                    ['Error when interpreting calibration file', ...
                    ' (LFR team''s RCT for BIAS/BICAS) "%s"'], ...
                    filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        
        
        
        % NOTE: TDS CWF cwfFactorsIvpt are already inverted (can be seen from
        % units).
        function RctData = read_TDS_CWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try                
                % NOTE: Undocumented in CDF: zVar CALIBRATION_TABLE is
                % volt/count for just multiplying the TDS signal (for this kind
                % of data). Is not a frequency-dependent transfer function.
                
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely on dataobj's
                % special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                
                factorsIvpt = shiftdim(Do.data.CALIBRATION_TABLE.data);
                
                % ASSERTIONS: Check CDF array sizes, no change in format.
                irf.assert.sizes(factorsIvpt, [3,1])
                
                RctData = [];
                RctData.factorsIvpt = factorsIvpt;
                
            catch Exc1
                Exc2 = MException(...
                    'BICAS:FailedToReadInterpretRCT', ...
                    ['Error when interpreting calibration file (TDS team''s', ...
                    ' LFM CWF RCT for BIAS/BICAS) "%s"'], ...
                    filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        
        
        
        % NOTE: The TDS RSWF RCT contains ITFs, not FTFs.
        %
        function RctData = read_TDS_RSWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special
                % behaviour for 1 record case.
                % ==> Remove leading singleton dimensions (just to be
                % format-tolerant), many assertions.
                freqsHz  = shiftdim(Do.data.CALIBRATION_FREQUENCY.data);   % 1x512 --> 512x1
                amplIvpt = shiftdim(Do.data.CALIBRATION_AMPLITUDE.data);   % 3x512 --> 3x512
                phaseDeg = shiftdim(Do.data.CALIBRATION_PHASE.data);       % 3x512 --> 3x512
                
                % ASSERTIONS: Check CDF array sizes, no change in format.
                nFreqs = irf.assert.sizes(...
                    freqsHz,  [-1,  1], ...
                    amplIvpt, [ 3, -1], ...
                    phaseDeg, [ 3, -1]);
                assert(nFreqs >= bicas.RCT.TF_TABLE_MIN_LENGTH)

                for iBlts = 1:3
                    % NOTE: RCT contains ITF, not FTF.
                    ItfIvpt = irf.utils.tabulated_transform(...
                        freqsHz * 2*pi, ...
                        amplIvpt(        iBlts, :), ...
                        deg2rad(phaseDeg(iBlts, :)));

                    % ASSERTION: INVERTED TF
                    assert(~ItfIvpt.toward_zero_at_high_freq(), ...
                        ['TDS RSWF transfer function appears to go toward', ...
                        ' zero at high frequencies. Has it not been', ...
                        ' inverted/made backward in time, i.e. does it not', ...
                        ' describe physical output-to-physical input?'])

                    ItfIvptList{iBlts} = ItfIvpt;
                end

                RctData = [];
                RctData.ItfIvptList = ItfIvptList;
                
            catch Exc1
                Exc2 = MException(...
                    'BICAS:FailedToReadInterpretRCT', ...
                    'Error when interpreting calibration file (TDS team''s', ...
                        ' LFM RSWF RCT for BIAS/BICAS) "%s"', ...
                    filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end



    end    % methods(Static, Access=public)
    
    
    
    methods(Static, Access=private)



        % Utility function
        %
        % Function for normalizing the indices of dataobj zVariables. dataobj
        % zVariable arrays have different meanings for their indices depending
        % on whether there are one record or many. If there is one record, then
        % there is not record index. If there are multiple records, then the
        % first index represents the record number. This function inserts a
        % size-one index as the first index.
        % 
        % DO   = dataobj(...)
        % data = Do.data.TRANSFER_FUNCTION_COEFFS.data
        %
        % NOTE: Not well tested on different types of zvar array sizes.
        % 
        function data = normalize_dataobj_ZV(DataobjZv)
            % PROPOSAL: Move to utils.
            % PROPOSAL: Shorter name:
            %   norm_dataobj_zv
            %   norm_do_zv
            
            data = DataobjZv.data;
            
            if DataobjZv.nrec == 1
                data = shiftdim(data, -1);
            end
        end



        % Utility function to simplify read_BIAS_RCT. Arguments correspond to
        % zVariables in BIAS RCT.
        %
        %
        % ARGUMENTS
        % =========
        % ftfNumCoeffs,
        % ftfDenomCoeffs : 2D matrix of numerator/denominator coefficients for
        %                  a sequence of FTFs. (iTime, iCoeff).
        %
        %
        % RETURN VALUE
        % ============
        % FtfArray : 1D column cell array of FTFs
        %            (irf.utils.rational_func_transform).
        %
        function FtfArray = create_TF_sequence(ftfNumCoeffs, ftfDenomCoeffs)
            
            % ASSERTIONS
            nTime = irf.assert.sizes(...
                ftfNumCoeffs,   [-1, -2], ...
                ftfDenomCoeffs, [-1, -2]);
            %assert(size(ftfNumCoeffs, 1) == size(ftfDenomCoeffs, 1))
            % The last FTF denominator coefficient (highest index, for which the
            % value is non-zero) must be =1.
            assert(...
                ftfDenomCoeffs(find(ftfDenomCoeffs, 1, 'last')) == 1, ...
                'BICAS:FailedToReadInterpretRCT', ...
                ['RCT should contain forward transfer function (FTF)', ...
                ' denominator coefficients,', ...
                ' where the highest-order (non-zero) coefficient', ...
                ' is the number one (1).', ...
                ' The data does not satisfy this criterion.'])

            FtfArray = {};
            for iTime = 1:nTime
                
                Ftf = irf.utils.rational_func_transform(...
                    ftfNumCoeffs(  iTime, :), ...
                    ftfDenomCoeffs(iTime, :));
                
                % ASSERTIONS
                assert(Ftf.has_real_impulse_response())
                % Assert FTF. Can not set proper error message.
                assert(Ftf.zero_in_high_freq_limit(), ...
                    'BICAS:FailedToReadInterpretRCT', ...
                    ['Transfer function is expected to be "forward",', ...
                    ' i.e. in the direction of the physical signal.', ...
                    ' It seems not to be.'])
                
                FtfArray{end+1, 1} = Ftf;    % Force column array.
            end
        end

        

    end    % methods(Static, Access=private)

    
    
end
