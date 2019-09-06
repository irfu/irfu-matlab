classdef calib
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR, TDS-CWF or TDS-RSWF) is identical (in all
% relevant parts) for RODP and ROC-SGSE.
%
%
% DEFINITIONS
% ===========
% CPV = counts/volt
% VPC = volt/count
% Rps = Radians/second
%
%
% UNFINISHED. NOT USED BY MAIN PROGRAM YET. NOT ENTIRELY CLEAR WHAT IT SHOULD INCLUDE.
%
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-15


% BOGIQ:
% ------
% PROPOSAL: Need multiple functions
%   Apply TF in frequency domain, in time domain(?)   (in both directions = inverse & non-inverse)
%   Convert TF from frequency domain to time domain (and maybe reverse for debugging)
%   Invert TF function.
%   Multiply/combine transfer functions.
%   Convert files data to time-domain transfer functions?
%   (Select one BIAS TF of many? Interpolate between transfer functions?!!!)
%   Select a TDS/LFR TF
%   Combine TFs: TDS/LFR, BIAS, possibly capacitance-TF.
%   Apply TF (inverse)
%   Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
%   Special treatment of 0 Hz?!
%
% QUESTION: How implement the parasitic capacitance TF? Use analytical TF directly, or convert to TF table?
% TODO: Read RCT
%
% PROPOSAL: apply_transfer_function_in_freq as separate utility function?
% PROPOSAL: Change name.
%   PRO: "calibration" is too long.
%   PRO: "calibration" overlaps with function names.
%   PROPOSAL: "calib", "cal".
%
% PROPOSAL: Function for converting frequency list, amplitude list, phase list, into struct (omega, Z).
%
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Function for permuting indices to handle dataobj's handling of 1 record-case.
%
% BOGIQ: RCT-reading functions
% ============================
% PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
% PROPOSAL: Assert CDF skeleton/master version number.
% PROPOSAL: Assert skeleton/master.
% PROPOSAL: Convert from file units to "mathematical units"? (radians/s, radians (for phase shift))
% PROPOSAL: Convert file transfer function to the right direction (invert it).
%   NOTE: Currently uncertain if such case exists, but it might.
% PROPOSAL: Assert/warn (depending on setting?) file units.
% PROPOSAL: Only use units in variable names.
% PROPOSAL: Use utility function for reading every zVariable.
%   PROPOSAL: Assert units from zVar attributes.



    properties(Access=private)
        
        % Different TFs. Initialized in constructor.
        % IMPLEMENTATION NOTE: Variable names include the prefix "BIAS" to indicate that these are the TFs for BIAS
        % as opposed to TDS and LFR which will have to be added later.
        %biasDcAbsGain        % BIAS gain for DC absolute (non-diff) data (instead of TF).
        %biasDcDiffGain       % BIAS gain for DC diff data (instead of TF).
        %biasDcAbsOffsets     % Unclear definition so far.
        
        %biasAcLoGainTf
        %biasAcHiGainTf
        
        %BiasDcRecordList = [];
        
        LfrTfTableList;
        tdsCwfFactors;
        TdsRswfTfTableList;
    end

    %###################################################################################################################

    methods(Access=public)

        function obj = calib(calibrationDir, pipelineId)
            % TODO-DECISION: Is it wise to specify the paths in the constructor? Read the filenames (and relative directory) from the constants instead?
            
            %obj.biasDcAbsGain  = 
            %obj.biasDcDiffGain = 
            
            %obj.BiasDcRecordList    = init_biasDcRecordList();
            
            filePath       = bicas.calib.find_RCT(calibrationDir, pipelineId, 'BIAS');
            bicas.calib.read_BIAS_RCT(filePath);
            
            filePath           = bicas.calib.find_RCT(calibrationDir, pipelineId, 'LFR');
            obj.LfrTfTableList = bicas.calib.read_LFR_RCT(filePath);
            
            filePath          = bicas.calib.find_RCT(calibrationDir, pipelineId, 'TDS-CWF');
            obj.tdsCwfFactors = bicas.calib.read_TDS_CWF_RCT(filePath);
            
            filePath               = bicas.calib.find_RCT(calibrationDir, pipelineId, 'TDS-RSWF');
            obj.TdsRswfTfTableList = bicas.calib.read_TDS_RSWF_RCT(filePath);
            
        end
        

        
%         function y2 = calibrate_LFR_DC_single(obj, dt, y1)
%             y2 = bicas.utils.apply_transfer_function(dt, y1, tfOmega, tfZ);
%             %y2 = y2 + 
%         end
% 
%         function calibrate_LFR_DC_diff(obj, )
%         end
% 
%         function calibrate_LFR_AC_low_gain(obj, )
%         end
% 
%         function calibrate_LFR_AC_high_gain(obj, )
%         end

    end    % methods(Access=public)
    
    %###################################################################################################################

    methods(Static, Access=public)
        
        
        
        % Evaluate (complex) analytical transfer function for arbitrary number of omega values.
        % 
        %     nc(1)*s^0 + ... + nc(nNc)*s^(nNc-1)
        % Z = -----------------------------------,   s = i*omega
        %     dc(1)*s^0 + ... + dc(nDc)*s^(nDc-1)
        %
        % nc = numeratorCoeffs
        % dc = denominatorCoeffs
        %
        % NOTE: MATLAB's "polyval" uses the opposite convention for the order/index of coefficients.
        %
        function Z = eval_analytical_transfer_func(omega, numeratorCoeffs, denominatorCoeffs)
            % PROPOSAL: Reclassify as generic function and move it.
            
            % ASSERTIONS: Required for using "flipud".
            assert(iscolumn(numeratorCoeffs))
            assert(iscolumn(denominatorCoeffs))
            % ASSERTION: Try to detect accidently switching the coefficients.
            % Denominator polynomial should have at least as high a degree as the numerator polynomial. <==> Z should
            % not diverge as omega-->inf.
            nNc = find(numeratorCoeffs,   1, 'last' );
            nDc = find(denominatorCoeffs, 1, 'last' );
            assert(nDc >= nNc)   
            
            s = 1i * omega;
            Z = polyval(flipud(numeratorCoeffs), s) ./ polyval(flipud(denominatorCoeffs), s);
        end
        
        
        
        function [BiasTfs, BiasCurrent, BiasVOffsets, BiasEOffsets] = read_BIAS_RCT(filePath)
            % NOTE: UNFINISHED. DOES NOT SAVE/RETURN OFFSETS.
            
            % TODO-DECISION: How handle time?
            %   PROPOSAL: "Only" access the BIAS values (trans.func and other) through a function instead of selecting indices in a data struct.
            %       PROPOSAL: (private method) [omegaRps, zVpc] = get_transfer_func(epoch, signalType)
            %           signalType = 'DC single' etc
            
            Do = dataobj(filePath);
            
            % Constants for interpreting the array indices in the CDF.
            NUMERATOR   = 1;
            DENOMINATOR = 2;
            DC_SINGLE = 1;
            DC_DIFF   = 2;
            AC_LG     = 3;
            AC_HG     = 4;
            
            try
                % NOTE: Assumes 1 CDF record or many (time-dependent values).
                % ==> Must handle that dataobj assigns differently for these two cases.
                Epoch_L                  = bicas.calib.norm_do_zv(Do.data.Epoch_L);
                Epoch_H                  = bicas.calib.norm_do_zv(Do.data.Epoch_H);
                BIAS_CURRENT_OFFSET      = bicas.calib.norm_do_zv(Do.data.BIAS_CURRENT_OFFSET);      % DEPEND_0 = Epoch_L
                BIAS_CURRENT_GAIN        = bicas.calib.norm_do_zv(Do.data.BIAS_CURRENT_GAIN);        % DEPEND_0 = Epoch_L
                V_OFFSET                 = bicas.calib.norm_do_zv(Do.data.V_OFFSET);                 % DEPEND_0 = Epoch_H
                E_OFFSET                 = bicas.calib.norm_do_zv(Do.data.E_OFFSET);                 % DEPEND_0 = Epoch_H
                TRANSFER_FUNCTION_COEFFS = bicas.calib.norm_do_zv(Do.data.TRANSFER_FUNCTION_COEFFS); % DEPEND_0 = Epoch_L
                
                nEpochL = size(Epoch_L, 1);
                nEpochH = size(Epoch_H, 1);
                
                % 1 CDF record : cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       F/TTT"   # 3=number of dimensions.
                % 2 CDF records: cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       T/TTT"

                % 1 CDF record:   size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [  4 2 8]
                % 2 CDF records:  size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [2 4 2 8]
                
                % size(TRANSFER_FUNCTION_COEFFS)              == [1  4  2  8]
                TRANSFER_FUNCTION_COEFFS = permute(TRANSFER_FUNCTION_COEFFS, [1, 4,3,2]);
                assert(size(TRANSFER_FUNCTION_COEFFS, 1) == nEpochL)
                assert(size(TRANSFER_FUNCTION_COEFFS, 2) >= 8)
                assert(size(TRANSFER_FUNCTION_COEFFS, 3) == 2)
                assert(size(TRANSFER_FUNCTION_COEFFS, 4) == 4)
                
                

                BiasTfs.dcSingle.numeratorCoeffs     = TRANSFER_FUNCTION_COEFFS(:, :, NUMERATOR,   DC_SINGLE);
                BiasTfs.dcSingle.denominatorCoeffs   = TRANSFER_FUNCTION_COEFFS(:, :, DENOMINATOR, DC_SINGLE);
                
                BiasTfs.dcDiff.numeratorCoeffs       = TRANSFER_FUNCTION_COEFFS(:, :, NUMERATOR,   DC_DIFF);
                BiasTfs.dcDiff.denominatorCoeffs     = TRANSFER_FUNCTION_COEFFS(:, :, DENOMINATOR, DC_DIFF);
                
                BiasTfs.acLowGain.numeratorCoeffs    = TRANSFER_FUNCTION_COEFFS(:, :, NUMERATOR,   AC_LG);
                BiasTfs.acLowGain.denominatorCoeffs  = TRANSFER_FUNCTION_COEFFS(:, :, DENOMINATOR, AC_LG);
                
                BiasTfs.acHighGain.numeratorCoeffs   = TRANSFER_FUNCTION_COEFFS(:, :, NUMERATOR,   AC_HG);
                BiasTfs.acHighGain.denominatorCoeffs = TRANSFER_FUNCTION_COEFFS(:, :, DENOMINATOR, AC_HG);
                
                BiasCurrent.offsets = BIAS_CURRENT_OFFSET;
                BiasCurrent.gains   = BIAS_CURRENT_GAIN;
                BiasVOffsets        = V_OFFSET;
                BiasEOffsets.E12 = E_OFFSET(:,1);
                BiasEOffsets.E13 = E_OFFSET(:,2);
                BiasEOffsets.E23 = E_OFFSET(:,3);
                
            catch Exc
                error('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath)
            end
        end



        % LfrTfTableList : LfrTfTableList{iFreq}{iBiasChannel}, iFreq=1..4 for F0..F3,
        %                   iFreq=1..3 : iBiasChannel=1..5 for BIAS_1..BIAS_5
        %                   iFreq=4    : iBiasChannel=1..3 for BIAS_1..BIAS_3
        %                  NOTE: This is different from LFR zVar FREQ.
        function LfrTfTableList = read_LFR_RCT(filePath)
            % BUG: 5+5+5+3 TFs, NOT 1+1+1+1.
            
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                
                freqsHz{1}  = shiftdim(Do.data.Freqs_F0.data);    % NOTE: Index {iLfrFreq}.
                freqsHz{2}  = shiftdim(Do.data.Freqs_F1.data);
                freqsHz{3}  = shiftdim(Do.data.Freqs_F2.data);
                freqsHz{4}  = shiftdim(Do.data.Freqs_F3.data);
                
                amplCpv{1}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F0.data);
                amplCpv{2}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F1.data);
                amplCpv{3}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F2.data);
                amplCpv{4}  = shiftdim(Do.data.TF_BIAS_123_amplitude_F3.data);
                
                phaseDeg{1} = shiftdim(Do.data.TF_BIAS_12345_phase_F0.data);
                phaseDeg{2} = shiftdim(Do.data.TF_BIAS_12345_phase_F1.data);
                phaseDeg{3} = shiftdim(Do.data.TF_BIAS_12345_phase_F2.data);
                phaseDeg{4} = shiftdim(Do.data.TF_BIAS_123_phase_F3.data);
                
                for iLfrFreq = 1:4
                    if iLfrFreq ~= 4
                        nBiasChannels = 5;
                    else
                        nBiasChannels = 3;
                    end

                    % ASSERTIONS: Check CDF array sizes, that no change in format.
                    assert(iscolumn(freqsHz{iLfrFreq}))
                    assert(ndims(amplCpv)  == 2)
                    assert(ndims(phaseDeg) == 2)
                    assert(size( amplCpv{iLfrFreq}, 2) == nBiasChannels)
                    assert(size(phaseDeg{iLfrFreq}, 2) == nBiasChannels)

                    for iBiasChannel = 1:nBiasChannels
                        omegaRps = freqsHz{iLfrFreq} * 2*pi;    % NOTE: Duplicates the same frequency list.
                        zVpc     = amplCpv{iLfrFreq}(:,iBiasChannel) .* exp(1i * deg2rad(phaseDeg{iLfrFreq}(:,iBiasChannel)));
                        
                        assert(isvector(omegaRps))
                        assert(isvector(zVpc))
                        assert(numel(omegaRps) == numel(zVpc))
                        
                        LfrTfTableList{iLfrFreq}{iBiasChannel}.omegaRps = omegaRps;
                        LfrTfTableList{iLfrFreq}{iBiasChannel}.zVpc     = zVpc;
                    end
                end
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        
        
        
        function tdsCwfFactorsVpc = read_TDS_CWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try                
                % NOTE: Undocumented in CDF: zVar CALIBRATION_TABLE is volt/count for just multiplying the TDS signal (for
                % this kind of data). Is not a frequency-dependent transfer function.
                
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                
                tdsCwfFactorsVpc = shiftdim(Do.data.CALIBRATION_TABLE.data);
                
                % ASSERTIONS: Check CDF array sizes, that no change in format.
                assert(iscolumn(tdsCwfFactorsVpc))
                assert(size(    tdsCwfFactorsVpc, 1) == 3)
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end            
        end
        
        
        
        function TdsRswfTfTableList = read_TDS_RSWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                freqsHz  = shiftdim(Do.data.CALIBRATION_FREQUENCY.data);
                amplVpc  = shiftdim(Do.data.CALIBRATION_AMPLITUDE.data);
                phaseDeg = shiftdim(Do.data.CALIBRATION_PHASE.data);
                
                % ASSERTIONS: Check CDF array sizes, that no change in format.
                assert(iscolumn(freqsHz));                
                assert(ndims(amplVpc)     == 2)
                assert(size( amplVpc, 1)  == 3)
                assert(ndims(phaseDeg)    == 2)
                assert(size( phaseDeg, 1) == 3)
                
                EJ_library.utils.assert.all_equal([...
                    length(freqsHz), ...
                    size(amplVpc,  2), ...
                    size(phaseDeg, 2) ]);
                
                for iBiasChannel = 1:3
                    omegaRps = freqsHz;    % NOTE: Duplicates the same frequency list.
                    zVpc     = shiftdim(amplVpc(iBiasChannel, :) .* exp(1i * deg2rad(phaseDeg(iBiasChannel, :))));
                    
                    assert(iscolumn(omegaRps))
                    assert(iscolumn(zVpc))
                    
                    TdsRswfTfTableList(iBiasChannel).omegaRps = omegaRps;
                    TdsRswfTfTableList(iBiasChannel).zVpc     = zVpc;
                end
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        

        
        % Find the path to the RCT to use, using the filenames specified in the documentation. If there are multiple
        % matching candidates, choose the latest one as indicated by the filename.
        %
        % RCT filenaming convention is described in:	ROC-PRO-DAT-NTT-00006-LES, 1/1 draft, Sect 4.3.2-3.
        %
        %
        % IMPLEMENTATION NOTE: BICAS is not following the documented filenaming convention here since (1) TDS and
        % LFR teams do not seem to use it and (2) it is uncertain how it can be applied to BIAS RCTs (which receiver
        % should the BIAS RCT specify?). Therefore using a regular expression that covers more files.
        %     NOTE: LFR RCTs use 2+12 digits instead of 10 in the timestamps (they add seconds=2 digits).
        %     NOTE: TDS RCTs use 2+6  digits instead of 10 in the timestamps (the have no time of day, only date)
        %
        % IMPLEMENTATION NOTE: Useful to have this as separate functionality so that the chosen RCT to use can be
        % explicitly overridden via e.g. settings.
        %
        %
        % ARGUMENTS
        % =========
        % pipelineId, rctId : String constants representing pipeline and RCT to be read.
        %
        function path = find_RCT(calibrationDir, pipelineId, rctId)            
            % PROPOSAL: Different regexp (not just the middle substring) for different types of RCT to handle different
            % naming conventions (e.g. different number of digits, whether to include substrings RCT, RPW).
            % PROPOSAL: Put type-dependent regexp. naming conventions in settings.
            %   PROPOSAL: Exclude pipelineId prefix and CDF suffix.
            
            % Examples of RCT filenames
            % -------------------------
            % BIAS:
            %       ROC-SGSE_CAL_RCT-BIAS_V201803211625.cdf   (old implemented convention)
            %       ROC-SGSE_CAL_RPW_BIAS_V201908231028.cdf   (new implemented convention, closer to documentation)
            %           SOLO_CAL_RCT-BIAS_V201901141146.cdf   (old implemented convention)
            % LFR:
            %       ROC-SGSE_CAL_RCT-LFR-BIAS_V20180724165443.cdf
            %           SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
            % TDS:
            %           SOLO_CAL_RCT-TDS-LFM-CWF-E_V20190128.cdf
            %           SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20190128.cdf
            %         (Two types of calibration files, but only RODP versions)

            %============================
            % Create regexp for filename
            %============================
            switch(pipelineId)
                case {'ROC-SGSE', 'RGTS'}
                    filenamePrefix = 'ROC-SGSE';
                case 'RODP'
                    filenamePrefix = 'SOLO';
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal pipelineId="%s"', pipelineId)
            end
            clear pipelineId
            
            % IMPLEMENTATION NOTE: This switch statement
            % (1) verifies the input/output, AND
            % (2) concentrates BICAS' knowledge of RCT filenaming to this function (the caller does not need to know
            %     anything about filenames).
            switch(rctId)
                case 'BIAS'
                    %rctDependentSubRegexp = 'RCT-BIAS_V20[0-9]{10}';   % Old convention.
                    rctDependentSubRegexp = 'RPW_BIAS_V20[0-9]{10}';
                case 'LFR'
                    rctDependentSubRegexp = 'RCT-LFR-BIAS_V20[0-9]{12}';
                case 'TDS-CWF'
                    rctDependentSubRegexp = 'RCT-TDS-LFM-CWF-E_V20[0-9]{6}';
                case 'TDS-RSWF'
                    rctDependentSubRegexp = 'RCT-TDS-LFM-RSWF-E_V20[0-9]{6}';
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId)
            end
            clear rctId
            
            %filenameRegexp = sprintf('%s_CAL_[A-Z_-]*%s_V20[0-9]{6,12}.(cdf|CDF)', filenamePrefix, rctDependentSubRegexp);
            filenameRegexp = sprintf('%s_CAL_%s.(cdf|CDF)', filenamePrefix, rctDependentSubRegexp);



            %=================================================
            % Find candidate files and select the correct one
            %=================================================
            dirObjectList = dir(calibrationDir);
            dirObjectList([dirObjectList.isdir]) = [];    % Eliminate directories.
            filenameList = {dirObjectList.name};
            filenameList(~EJ_library.utils.regexpf(filenameList, filenameRegexp)) = [];    % Eliminate non-matching filenames.
            
            % ASSERTION / WARNING
            if numel(filenameList) == 0
                % ERROR
                error('BICAS:calib:CannotFindRegexMatchingRCT', ...
                    'Can not find any calibration file that matches regular expression "%s" in directory "%s".', ...
                    filenameRegexp, calibrationDir);
            end
            % CASE: There is at least on candidate file.
            
            filenameList = sort(filenameList);
            filename     = filenameList{end};
            path = fullfile(calibrationDir, filename);
            
            if numel(filenameList) > 1
                % WARNING/INFO/NOTICE
                msg = sprintf('Found multiple calibration files matching regular expression "%s"\nin directory "%s".\nSelecting the latest one as indicated by the filename: "%s".\n', ...
                    filenameRegexp, calibrationDir, filename);
%                 for i = 1:numel(filenameList)
%                     msg = [msg, sprintf('    %s\n', filenameList{i})];
%                 end
                bicas.log('debug', msg)
            end
            
            % IMPLEMENTATION NOTE: Not logging which calibration file is selected, since this function is not supposed
            % to actually load the content.
        end



%         function tfZ = parasitic_capacitance_TF(tfOmega)
%             % Calculate Z(omega) values for TF representing parasitic capacitances (based on analytic function).
% 
%             % Function name? "Input capacitance"?
%             % Not read R & C from constants here? Submit as arguments?
%             capacitanceFarad =
%             impedanceOhm     =
%             
%             % Correct for a TF?!
%             tfZ = 1 / (1 + 1j*tfOmega*capacitanceFarad*impedanceOhm);
%             
%             error('BICAS:calib:OperationNotImplemented', 'Function not implemented Yet.')
%         end



%         function biasCurrentAmpere = HK_TM_bias_to_bias_current(hkBiasTm)
%         % Convert HK diagnostic TM bias current values to physical units.
%         %
%         % NOTE: The HK bias current values are the measured onboard analog bias current values but are only meant as
%         % diagnostic values, not the proper bias current values. Therefore the calibrated values should only be seen as
%         % approximate.
%         % NOTE: The conversion function can be found in the BIAS spec, under BIAS1/2/3 under "Telemetry". (Not to be confused with the corresponding telecommands.) 
%         % NOTE: The conversion functions are identical for all three probes.
%
%         % NOTE: 
%         end



    end    % methods(Static, Access=public)
    
    
    
    methods(Static, Access=private)
        
        
        
        % Function for normalizing the indices of dataobj zVariables.
        % dataobj zVariable arrays have different meanings for their indices depending on whether there are one record
        % or many. If there is one record, then there is not record index. If there are multiple records, then the first
        % index represents the record number. This function inserts a size-one index as the first index.
        % 
        % Do = dataobj(...)
        % data = Do.data.TRANSFER_FUNCTION_COEFFS.data
        %
        % NOTE: Not well tested on different types of zvar array sizes.
        % 
        function data = norm_do_zv(DataobjZVar)
            % PROPOSAL: Move to utils.
            % PROPOSAL: Shorter name:
            %   norm_dataobj_zvar
            %   norm_do_zv_data
            %   norm_do_zv
            
            data = DataobjZVar.data;            
            
            if DataobjZVar.nrec == 1
                %nDims = ndims(data);
                %order = [nDims + 1, 1:nDims];
                %data = permute(data, order);
                data = shiftdim(data, -1);
            end
        end

            

    end    %methods(Static, Access=private)

end


        
        % TEST
%         function inputSignalVolt = calibrate_DC(obj, temperatureCelsius, stimuli, antennaCh, lfrCh, outputSignalVolt)
%             %if numel(antennaCh) == 2
%             %    antennaCh = sort(antennaCh);
%             %end
%             
%             if     isequal(antennaCh, [1])
%             elseif isequal(antennaCh, [2])
%             elseif isequal(antennaCh, [3])
%             elseif isequal(antennaCh, [1,2])
%             elseif isequal(antennaCh, [1,3])
%             elseif isequal(antennaCh, [1,3])
%             else
%                 error('BICAS:Assertion', 'Can not interpret antenna channels.')
%             end
%             
%             inputSignalVolt = offset + slope .* outputSignalVolt;
%         end
        