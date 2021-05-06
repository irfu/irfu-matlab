%
% Class that collects "processing functions" as public static methods. Only
% covers processing L2-->L3 ("23" in the name).
%
% This class is not meant to be instantiated.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% See bicas.calib.
% ZV   : CDF zVariable, or something analogous to it. If refers to CDF:ish
%        content, then the first index corresponds to the CDF record.
% SPR  : Samples Per (CDF) Record. Only refers to actual data (currents,
%        voltages), not metadata.
% UFV  : Use Fill Values (refers to records which data should overwritten with
%        fill values)
% ORIS : Oiginal sampling (frequency), as opposed to DWNS.
% DWNS : Downsampled, as opposed to ORIS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef proc_sub23
%##############################################################################################
%
% PROPOSAL: POLICY: Include all functions which set "policy"/configure the output of datasets.
%
% PROPOSAL: Test code for code that downsamples.
%   Ex: bicas.proc_sub23.downsample_bin_L12_QUALITY_BITMASK() -- Too trivial?
%   Ex: bicas.proc_sub23.downsample_bin_QUALITY_FLAG()        -- Too trivial?
%   Ex: bicas.proc_sub23.downsample_bin_sci_values()          -- Already has test code
%   Ex: bicas.proc_sub23.downsample_Epoch()                   -- Already has test code
%   --
%   PRO: Can verify now uncertain edge cases.
%       Ex: Quality zVars when science data=fill values.
%       PRO: Can verify future bugfix for integer quality zVar=fill value when
%            there is no science data.
%
% PROPOSAL: Split up processing between (a) density, and (b) E-field & SCPOT.
%   PRO: Faster
%       CON: Not very heavy operation.
%   PRO: Leads to better organization of code.
%       PRO: process_L2_to_L3() is too large and should be split up anyway.
%
% PROPOSAL: QUALITY_FLAG and L2_QUALITY_BITMASK should be FILL VALUE (not zero) for empty
%           bins.
%   PROBLEM: BICAS might not support that behaviour yet. /2021-05-06
%
%##############################################################################################



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################

    
    
        % Processing function for processing L2-->L3 (not VHT).
        %
        function [OutEfieldOris,  OutEfieldDwns, ...
                  OutScpotOris,   OutScpotDwns, ...
                  OutDensityOris, OutDensityDwns] ...
                = process_L2_to_L3(InLfrCwf, SETTINGS, L)

            % PROPOSAL: Split up in one part for non-downsampled and
            %           downsampled.
            %
            % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
            %           (still combine non-downsampled and downsampled).
            %   CON: Slows down overall processing.
            %       PRO: Must read same L2 dataset multiple times.
            %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
            %   CON: There is much shared functionality for 3 quality zVars.
            %       PRO: Same ~constants
            %           Ex: INPUT_DATASET_ID, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
            %       PRO: Read setting QUALITY_FLAG_MIN_FOR_USE
            %       PRO: Normalizing CWF zVar names.
            %       PRO: Preparations for downsampled.
            %           Bin locations, bundling of records,
            %           Downsampling of quality variables
            %               (QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK).
            %           DELTA_PLUS_MINUS_dwns
            %
            % BUG: Fill values in the __INPUT__
            %   QUALITY_FLAG,
            %   QUALITY_BITMASK,
            %   L2_QUALITY_BITMASK are not known/recognized since the variables
            %   are not double.
            %   NOTE: Unrelated ROC BUG:
            %               https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/48
            %         L1 QUALITY_BITMASK seems to use the wrong value (255) as
            %         fill value (FILLVAL=65535). ==> A bug fix would not fix
            %         the entire issue.
            %   PROPOSAL: Use double also for CDF integer variables so NaN can
            %             represent fill value also for these.
            %
            % NOTE: L2 LFR-CWF-E skt had zVar
            %   QUALITY_BITMASK=CDF_UINT1, fill value=255 (wrong)
            % until skt V12 when it was changed to
            %   QUALITY_BITMASK=CDF_UINT2, fill value 65535 (correct).


            
            %===========
            % Constants
            %===========
            % The only acceptable input DATASET_ID.
            INPUT_DATASET_ID          = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            BIN_LENGTH_WOLS_NS        = int64(10e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);
            % Regular expression for the format of version strings from
            % BICAS-external code.
            % Equivalent to: yyyy-mm-ddThh:mm:ss
            CODE_VER_STR_REGEXP = ...
                '[0-9]{4}-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]';
            


            %======================
            % Normalize zVar names
            %======================
            [InLfrCwf.Zv, fnChangeList] = ...
                EJ_library.utils.normalize_struct_fieldnames(InLfrCwf.Zv, ...
                {{{'VDC', 'V'}, 'VDC'}}, 'Assert one matching candidate');

            bicas.proc_utils.handle_zv_name_change(...
                fnChangeList, INPUT_DATASET_ID, SETTINGS, L, 'VDC', ...
                'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



            %=================================================================
            % Call BICAS-external code to calculate (EFIELD, SCPOT) + DENSITY
            %=================================================================
            R = bicas.proc_sub23.calc_EFIELD_SCPOT(InLfrCwf.Zv, SETTINGS);
            %
            [NeScpTs, psp2neCodeVerStr] = bicas.proc_sub23.calc_DENSITY(R.PspTs);



            %====================================================================
            % Derive values for CDF global attribute "Misc_calibration_versions"
            %====================================================================
            assert(isempty(R.vdccalMatVerStr), ...
                ['solo.vdccal() no longer returns an empty vdccalMatVerStr', ...
                ' (string representing the version of the corresponding', ...
                ' .mat file). BICAS therefore needs to be updated.'])
            EJ_library.assert.castring_regexp(R.vdccalCodeVerStr, CODE_VER_STR_REGEXP)
            EJ_library.assert.castring_regexp(psp2neCodeVerStr,   CODE_VER_STR_REGEXP)
            %
            % NOTE: Should not add BICAS version to glob.attr.
            % "Misc_calibration_versions" since this is already encoded in
            % global attribute "Software_version" (together with
            % "Software_name").
            %
            % NOTE: There version string for the solo.vdccal() .mat file has not
            % been implemented yet.
            % NOTE: Density Misc_calibration_versions contains both versions,
            % since density is derived from PSP.
            vdccalStr = ['solo.vdccal() code version ', R.vdccalCodeVerStr];
            psp2neStr = ['solo.psp2ne() code version ', psp2neCodeVerStr];
            gaEfieldScpot_Misc_calibration_versions = {vdccalStr};
            gaDensity_Misc_calibration_versions     = {vdccalStr, psp2neStr};



            %=========================================
            % Misc. variables shared between datasets
            %=========================================
            % Global attributes, shared between ALL 3x2 datasets.
            InitialGa = struct();
            InitialGa.OBS_ID                 = InLfrCwf.Ga.OBS_ID;
            InitialGa.SOOP_TYPE              = InLfrCwf.Ga.SOOP_TYPE;
            % zVariables, shared between all ORIS/NON-DOWNSAMPLED datasets.
            InitialOrisZv = struct();
            InitialOrisZv.Epoch              = InLfrCwf.Zv.Epoch;
            InitialOrisZv.QUALITY_BITMASK    = InLfrCwf.Zv.QUALITY_BITMASK;
            InitialOrisZv.L2_QUALITY_BITMASK = InLfrCwf.Zv.L2_QUALITY_BITMASK;
            InitialOrisZv.QUALITY_FLAG       = InLfrCwf.Zv.QUALITY_FLAG;
            InitialOrisZv.DELTA_PLUS_MINUS   = InLfrCwf.Zv.DELTA_PLUS_MINUS;
            %
            InitialOris = struct(...
                'Ga', InitialGa, ...
                'Zv', InitialOrisZv);
            clear InitialZv
            %
            [InitialDwnsZv, iRecordsInBinCa] = bicas.proc_sub23.init_shared_downsampled(...
                InLfrCwf.Zv, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS);
            InitialDwns = struct('Zv', InitialDwnsZv);
            % NOTE: Not setting DWNS .Ga here, since DWNS datasets later copy
            % .Ga from the respective ORIS datasets.
            


            %====================================
            % zVars for EFIELD (not downsampled)
            %====================================
            OutEfieldOris = InitialOris;
            OutEfieldOris.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutEfieldOris.Zv.EDC_SRF                   = R.zvEdcMvpm;



            %===================================
            % zVars for SCPOT (not downsampled)
            %===================================
            OutScpotOris = InitialOris;
            OutScpotOris.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutScpotOris.Zv.SCPOT                     = R.ScpotTs.data;
            OutScpotOris.Zv.PSP                       = R.PspTs.data;



            %=====================================
            % zVars for DENSITY (not downsampled)
            %=====================================
            OutDensityOris = InitialOris;
            OutDensityOris.Ga.Misc_calibration_versions = gaDensity_Misc_calibration_versions;
            %
            OutDensityOris.Zv.DENSITY                   = NeScpTs.data;



            %==============================
            % zVars for EFIELD DOWNSAMPLED
            %==============================
            OutEfieldDwns    = InitialDwns;
            OutEfieldDwns.Ga = OutEfieldOris.Ga;
            %
            [OutEfieldDwns.Zv.EDC_SRF, ...
             OutEfieldDwns.Zv.EDCSTD_SRF] = bicas.proc_sub23.downsample_sci_zVar(...
                OutEfieldOris.Zv.EDC_SRF, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa);
            %
            % NOTE: Merge across samples in same record.
            %b = all(isnan(OutEfieldDwns.Zv.EDC_SRF), 2);
            %OutEfieldDwns.Zv.QUALITY_FLAG(b) = 0;

            

            %=============================
            % zVars for SCPOT DOWNSAMPLED
            %=============================
            OutScpotDwns    = InitialDwns;
            OutScpotDwns.Ga = OutScpotOris.Ga;
            %
            [OutScpotDwns.Zv.SCPOT, ...
             OutScpotDwns.Zv.SCPOTSTD] = bicas.proc_sub23.downsample_sci_zVar(...
                OutScpotOris.Zv.SCPOT, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa);
            %
            [OutScpotDwns.Zv.PSP, ...
             OutScpotDwns.Zv.PSPSTD] = bicas.proc_sub23.downsample_sci_zVar(...
                OutScpotOris.Zv.PSP, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa);
            %
            %b = isnan(OutScpotDwns.Zv.SCPOT) & ...
            %    isnan(OutScpotDwns.Zv.PSP);
            %OutScpotDwns.Zv.QUALITY_FLAG(b) = 0;



            %===============================
            % zVars for DENSITY DOWNSAMPLED
            %===============================
            OutDensityDwns    = InitialDwns;
            OutDensityDwns.Ga = OutDensityOris.Ga;
            %
            [OutDensityDwns.Zv.DENSITY, ...
             OutDensityDwns.Zv.DENSITYSTD] = bicas.proc_sub23.downsample_sci_zVar(...
                OutDensityOris.Zv.DENSITY, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa);
            %
            %b = isnan(OutDensityDwns.Zv.DENSITY);
            %OutDensityDwns.Zv.QUALITY_FLAG(b) = 0;

        end    % process_L2_to_L3



        % Downsample a single NxM science zVar.
        %
        % Use bins and for every bin, derive median and modified standard
        % deviation over dimension 1 (within bin). Construct two zVariables for
        % median+MSTD for the corresponding downsampled CDF record.
        %
        % NOTE: Can handle zero records in bin.
        % NOTE: Function is only public so that automated test code can access
        % it.
        %
        %
        % ARGUMENTS
        % =========
        % zv
        %       (iCdfRecord, iChannel).
        % nMinReqSamples
        %       Minimum number of samples (fill value or not) for not returning
        %       fill value.
        % iRecordsInBinCa
        %       Column cell array. {iBin, 1}
        %
        %
        % RETURN VALUES
        % =============
        % zvMed          : (iBin, iChannel). Median
        % zvMstd         : (iBin, iChannel). Modified STandard Deviation (MSTD).
        % bTooFewRecords : (iBin, 1). Logical. Whether number of samples fell
        %                   below threshold for particular bin. Potentially
        %                   useful for setting QUALITY_FLAG. Not used(?).
        %
        function [zvMed, zvMstd, bTooFewRecords] = downsample_sci_zVar(...
                zv, nMinReqRecords, iRecordsInBinCa)
            
            % PROPOSAL: Incorporate bicas.proc_sub23.downsample_bin_sci_values()
            %           in this function.
            %   NOTE: Only used here.
            %   PRO: Faster?
            %   PRO: Clearer?
            %
            % PROPOSAL: Require nMinReqSamples >= 1? Code can handle 0, though it gives NaN.

            % ASSERTION
            assert(isfloat(zv))
            assert(nMinReqRecords >= 0)
            [nSpr, nRecordsDwns] = EJ_library.assert.sizes(...
                zv,              [NaN, -1], ...
                iRecordsInBinCa, [-2]);
            
            % Pre-allocate
            zvMed          = NaN(  nRecordsDwns, nSpr);
            zvMstd         = NaN(  nRecordsDwns, nSpr);
            bTooFewRecords = false(nRecordsDwns, 1);    % Default values.

            for iBin = 1:nRecordsDwns
                k           = iRecordsInBinCa{iBin};
                
                binZv       = zv(k, :);
                nBinRecords = size(binZv, 1);

                if nBinRecords < nMinReqRecords
                    % CASE: Too few samples.
                    binMed     = NaN;
                    binMstd    = NaN;                    
                    bTooFewRecords(iBin) = true;
                else

                    binMed  = median(binZv, 1);
                    binMstd = NaN(1, nSpr);    % Pre-allocate.
                    for j = 1:nSpr
                        binMstd(1, j) = bicas.utils.modif_std_deviation(...
                            binZv(:, j), ...
                            binMed(j), ...
                            1);
                    end
                end
                
                zvMed (iBin, :) = binMed;
                zvMstd(iBin, :) = binMstd;
                
                % clear binMed binMstd
            end
            
        end



        % NOTE: THIS FUNCTION HAS BEEN COPIED INTO downsample_sci_zVar().
        % OBSOLETE.
        % 
        % Derive median and modified standard deviation over dimension 1. For a
        % range of CDF records in a zVariable (at most 1D/CDF record), construct
        % two zVariables for median+MSTD for the corresponding downsampled CDF
        % record.
        %
        % NOTE: Can handle zero input records.
        % NOTE: Function is only public so that automated test code can access
        % it.
        %
        %
        % ARGUMENTS
        % =========
        % zvSegment
        %       (iCdfRecord, iChannel).
        % nMinReqSamples
        %       Minimum number of samples (fill value or not) for not returning
        %       fill value.
        %
        %
        % RETURN VALUES
        % =============
        % med  : (1, iChannel). 1xN. Median
        % mstd : (1, iChannel). 1xN. Modified STandard Deviation (MSTD).
        %
%         function [med, mstd] = downsample_bin_sci_values(...
%                 zvSegment, nMinReqSamples)
% 
%             % PROPOSAL: Argument for minimum number of samples in each bin. If
%             %           number of samples per bin is below limit, then return
%             %           NaN.
%             %   PROPOSAL: Take NaN samples into account. Exclude them.
%             %       CON: Can not do for
%             %
%             % PROPOSAL: Include the loop over downsampled records.
%             %   PRO: Same procedure for all downsampled datasets.
%             %   NOTE: Future variations in procedure could be handle using
%             %   parameters.
%             %   CON: ~Can/should still not eliminate setting quality zVariables
%             %        in loop.
%             %
%             % ASSERTION
%             % Only first two dimensions may be size non-one (with current
%             % implementation).
%             assert(ismatrix(zvSegment))
%             assert(isscalar(nMinReqSamples))
% 
% 
% 
%             nRecords = size(zvSegment, 1);
%             nSpr     = size(zvSegment, 2);   % SPR = Samples Per (CDF) Record
% 
%             % ~NORMALIZATION
%             if nRecords < nMinReqSamples
%                 % CASE: Too few samples. ==> Remove all samples. ==> MSTD=NaN
%                 zvSegment = zvSegment([], :);
%             end
% 
%             med  = median(zvSegment, 1);
%             mstd = NaN(1, nSpr);    % Pre-allocate.
%             for i = 1:nSpr
%                 mstd(1, i) = bicas.utils.modif_std_deviation(...
%                     zvSegment(:, i), med(i), 1);
%             end
%         end
        
        
        
        % Utility function for downsampling data by grouping together adjacent
        % time intervals that have the same length when discounting leap
        % seconds.
        %
        % NOTE: Function is only public so that automated test code can access
        % it.
        %
        %
        % ARGUMENTS
        % =========
        % zvAllTt2000           : Column array. ~Epoch.
        %                         PROBLEM: Can not handle zvAllTt2000(1),
        %                         zvAllTt2000(end) being during positive leap
        %                         second.
        % boundaryRefTt2000     : Must not be during leap second. Specifies
        %                         boundaries will be together with other
        %                         arguments.
        % binLengthWolsNs       : Length of each bin.
        % binTimestampPosWolsNs : Position of timestamp that represents bin,
        %                         relative to beginning of bin.
        %
        %
        % RETURN VALUES
        % =============
        % zvTt2000
        %       Column array. ~Epoch. One timestamp per bin.
        % iRecordsInBinCa
        %       Column cell array.
        %       Indices to CDF records for respective bins.
        %       {iBin, 1}(iSamples,1) = Non-downsampled CDF record number.
        % nRecordsPerBin
        %       Column array. (iBin, 1) = Number of non-downsampled records in
        %       bin. Could be useful for setting QUALITY_FLAG. Currently not
        %       used(?).
        % binSizeArrayNs
        %       (iBin, 1). Bin size.
        %       RATIONALE: Useful for automatic testing, setting zVar
        %       DELTA_PLUS_MINUS (if one wants to account for leap seconds).
        %
        %
        % NAMING CONVENTIONS
        % ==================
        % WOLS     = WithOut Leap Seconds
        % TTW      = TT2000 WOLS
        % bin      = Time interval within which all corresponding CDF records
        %             should be condensed to one.
        % boundary = Edge of bin(s).
        %
        function [zvBinsTt2000, iRecordsInBinCa, nRecordsPerBin, binSizeArrayNs] = ...
            downsample_Epoch(...
                zvAllTt2000, boundaryRefTt2000, ...
                binLengthWolsNs, binTimestampPosWolsNs)
            
            % PROPOSAL: Return boundariesTt2000 instead of binSizeArrayNs.
            %   PRO: More information.
            %   PRO: Easy to derive binSizeArrayNs = diff(boundariesTt2000);
            %   CON: Undefined (?) for special case zero bins.
            %
            % PROPOSAL: Argument for minimum number of samples in each bin. If
            %           number of samples per min is below limit, then they are
            %           excluded.
            %   CON: Samples may be NaN, but this function does not have access
            %        to that information.
            

            
            % ASSERTIONS
            bicas.proc_utils.assert_zv_Epoch(zvAllTt2000)
            % NOTE: Function algorithm assumes this monotonic increase.
            assert(issorted(zvAllTt2000, 'strictascend'))
            %
            bicas.proc_utils.assert_zv_Epoch(boundaryRefTt2000)
            assert(isscalar(boundaryRefTt2000))
            assert(isscalar(binLengthWolsNs))
            assert(isa(binLengthWolsNs,       'int64'))
            assert(isa(binTimestampPosWolsNs, 'int64'))
            assert((0 <= binTimestampPosWolsNs)...
                && (binTimestampPosWolsNs <= binLengthWolsNs))
            
            
            
            if isempty(zvAllTt2000)
                % CASE: zvAllTt2000 is empty.
                zvBinsTt2000    = int64(ones(0,1));
                iRecordsInBinCa = cell(0,1);
                nRecordsPerBin  = zeros(0,1);
                binSizeArrayNs  = zeros(0,1);
                return
            end
            % CASE: zvAllTt2000 is not empty.
            
            
            
            ttw1           = EJ_library.cdf.time.TT2000_to_TT2000WOLS(zvAllTt2000(1));
            ttw2           = EJ_library.cdf.time.TT2000_to_TT2000WOLS(zvAllTt2000(end));
            boundaryRefTtw = EJ_library.cdf.time.TT2000_to_TT2000WOLS(boundaryRefTt2000);
            
            %======================================
            % Find bin boundaries & bin timestamps
            %======================================
            % "Round" ttw1 down to nearest lower interval boundary.
            ttw1Floor = ...
                idivide(ttw1 - boundaryRefTtw, binLengthWolsNs, 'floor') ...
                * binLengthWolsNs + boundaryRefTtw;
            
            % Find smallest number of time intervals that will cover (and exceed
            % if necessary) ttw1 to ttw2.
            nBins         = idivide(ttw2 - ttw1Floor, binLengthWolsNs, 'ceil');
            
            boundariesTtw = (ttw1Floor + [0:nBins] * binLengthWolsNs)';
            zvBinsTtw     = boundariesTtw(1:end-1) + binTimestampPosWolsNs;
            
            
            boundariesTt2000 = EJ_library.cdf.time.TT2000WOLS_to_TT2000(boundariesTtw);
            zvBinsTt2000     = EJ_library.cdf.time.TT2000WOLS_to_TT2000(zvBinsTtw);
            binSizeArrayNs   = diff(boundariesTt2000);
            
            %==================
            % Assign iRecordCa
            %==================
            iRecordsInBinCa = cell( nBins, 1);
            nRecordsPerBin  = zeros(nBins, 1);
            for i = 1:nBins
                % Slow?
                % PROPOSAL: Speed up by using that zvTt2000 and boundaries are
                % sorted. Iterate over records, stopping at boundaries.
                iRecordsInBinCa{i} = find(...
                      (boundariesTt2000(i) <= zvAllTt2000)...
                    & (zvAllTt2000         <  boundariesTt2000(i+1)));
                
                nRecordsPerBin(i) = length(iRecordsInBinCa{i});
            end
        end
        


    end    % methods(Static, Access=public)



    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################
    
    
    
        % Calculate both
        %   (1) ELECTRIC FIELD, and
        %   (2) SPACECRAFT POTENTIALS
        % via the same BICAS-external code solo.vdccal() (still inside
        % irfu-matlab).
        %
        % Largely a wrapper around solo.vdccal().
        %
        % NOTE: Needs to be careful with the units, and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % RETURN VALUE
        % ============
        % R : Struct with multiple variables.
        %     NOTE: Return values are packaged as a struct to provide named
        %     return values and avoid confusing similar return results with each
        %     other.
        %
        function R = calc_EFIELD_SCPOT(...
                InLfrCwfZv, SETTINGS)
            
            QUALITY_FLAG_minForUse = SETTINGS.get_fv(...
                'PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN');


            
            % Shorten recurring variables.
            zv_VDC   = InLfrCwfZv.VDC;
            zv_EDC   = InLfrCwfZv.EDC;
            zv_Epoch = InLfrCwfZv.Epoch;            
            
            %======================================================
            % Create input variables for solo.vdccal()
            % ----------------------------------------
            % Set records to NaN for QUALITY_FLAG below threshold.
            %======================================================
            bDoNotUse = InLfrCwfZv.QUALITY_FLAG < QUALITY_FLAG_minForUse;
            zv_VDC(bDoNotUse, :) = NaN;
            zv_EDC(bDoNotUse, :) = NaN;
            %
            % NOTE: Should TSeries objects really use TensorOrder=1 and
            % repres={x,y,z}?!! VDC and EDC are not time series of vectors, but
            % fo three scalars. Probably does not matter. solo.vdccal() does
            % indeed use VDC.x, EDC.x etc.
            VdcTs = TSeries(...
                EpochTT(zv_Epoch), zv_VDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            EdcTs = TSeries(...
                EpochTT(zv_Epoch), zv_EDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            
            
            
            %==========================
            % CALL BICAS-EXTERNAL CODE
            %==========================
            [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] ...
                = solo.vdccal(VdcTs, EdcTs);
            clear VdcTs EdcTs
            %==========================
            
            
            
            % ASSERTIONS
            EJ_library.assert.sizes(...
                zv_Epoch,      [-1, 1], ...
                EdcSrfTs.data, [-1, 3], ...
                PspTs.data,    [-1, 1], ...
                ScpotTs.data,  [-1, 1]);
            assert(strcmp(EdcSrfTs.units,            'mV/m'))
            assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF'))
            assert(strcmp(PspTs.units,               'V'))
            assert(strcmp(ScpotTs.units,             'V'))
            
            
            
            %===================================================================
            % Normalize the representation of E-field X-component
            % (EdcSrfTs --> zvEdcMvpm)
            % ---------------------------------------------------
            % Set E_x = NaN, but ONLY if assertion deems that the corresponding
            % information is missing.
            %
            % IMPLEMENTATION NOTE: solo.vdccal set antenna 1 to be zero, if the
            % source data is non-fill value/NaN, but NaN if fill value. Must
            % therefore check for both zero and NaN.
            % Ex: Dataset 2020-08-01
            %===================================================================
            zvEdcMvpm = EdcSrfTs.data;    % MVPM = mV/m
            clear EdcSrfTs
            % IMPLEMENTATION NOTE: ismember() does not work for NaN.
            assert(all(zvEdcMvpm(:, 1) == 0 | isnan(zvEdcMvpm(:, 1))), ...
                ['EDC for antenna 1 returned from', ...
                ' solo.vdccal() is neither zero nor NaN and can therefore', ...
                ' not be assumed to be unknown anymore.', ...
                ' Verify that this is correct solo.vdccal() behaviour and', ...
                ' (if correct) then update BICAS to handle this.'])
            zvEdcMvpm(:, 1) = NaN;
            
            
            
            % Prepare return struct.
            R = [];
            R.PspTs            = PspTs;
            R.ScpotTs          = ScpotTs;
            R.zvEdcMvpm        = zvEdcMvpm;
            R.vdccalCodeVerStr = vdccalCodeVerStr;
            R.vdccalMatVerStr  = vdccalMatVerStr;
            
        end



        % Calculate DENSITY via a BICAS-external code solo.psp2ne() (still
        % inside irfu-matlab).
        %
        % Largely a wrapper around solo.psp2ne().
        % 
        % NOTE: Needs to be careful with the units, and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % NOTE: Empirically, some return values are NaN.
        % NOTE: Shortening "SCP" comes from the return variable name in
        % solo.psp2ne(). Do not know what it means.
        %
        function [NeScpTs, psp2neCodeVerStr] = calc_DENSITY(PspTs)
            
            %==========================
            % CALL BICAS-EXTERNAL CODE
            %==========================
            [NeScpTs, psp2neCodeVerStr] = solo.psp2ne(PspTs);
            %==========================
            
            % ASSERTIONS
            EJ_library.assert.sizes(...
                PspTs.data,   [-1, 1], ...
                NeScpTs.data, [-1, 1]);
            assert(all( (NeScpTs.data > 0) | isnan(NeScpTs.data)), ...
                'solo.psp2ne() returned non-positive (non-NaN) plasma density.')
            assert(strcmp(NeScpTs.units, 'cm^-3'))
            
        end



        % Derive values which are used by all DOWNSAMPLED datasets.
        %
        % RETURN VALUES
        % =============
        % InitialDwnsZv
        %       Struct with zVariables. 
        % iRecordsInBinCa
        %       Distribution of non-downsampled records in bins.
        %
        function [InitialDwnsZv, iRecordsInBinCa] = init_shared_downsampled(...
                InLfrCwfZv, binLengthWolsNs, binTimestampPosWolsNs)
            
            %================
            % Calculate bins
            %================
            % Find bin boundary reference timestamp. This is used for
            % setting the bin boundaries together with the bin length.
            v = spdfbreakdowntt2000(InLfrCwfZv.Epoch(1));
            % UTC subsecond (milliseconds, microseconds, nanoseconds)
            v(6)   = 5;   % UTC second
            v(7:9) = 0;
            boundaryRefTt2000 = spdfcomputett2000(v);
            % Find
            % (1) bin timestamps (downsampled timestamps to represent each bin),
            %     and
            % (2) which (non-downsampled) records belong to which bins
            %     (=downsampled records).
            [zvEpochDwns, iRecordsInBinCa, ~, binSizeArrayNs] = ...
                bicas.proc_sub23.downsample_Epoch(...
                    InLfrCwfZv.Epoch, ...
                    boundaryRefTt2000, ...
                    binLengthWolsNs, ...
                    binTimestampPosWolsNs);
            nRecordsDwns = numel(zvEpochDwns);
            
            
            
            %====================================================================
            % Derive downsampled versions of quality zVariables
            % -------------------------------------------------
            %
            % QUALITY_BITMASK and L2_QUALITY_BITMASK
            % --------------------------------------
            % QUALITY_BITMASK and L2_QUALITY_BITMASK should by their very
            % definition only be COPIED FROM LOWER ARCHIVING LEVELS (input
            % datasets). However, in practice they have to be modified somewhat
            % due to the downsampling, but this should under any circumstance be
            % independent of the data/data processing (non-quality variables)
            % itself.
            %
            % There are two cases which need to be handled/configured
            % intelligently:
            % (1) Calculate a combined value for every bin with AT LEAST ONE
            %     non-downsampled timestamp.
            % (2) Set value for every bin with ZERO non-downsampled timestamps.
            %     (Zero? fill value?)
            %
            %
            % QUALITY_FLAG
            % ------------
            % IMPORTANT: QUALITY_FLAG needs to be downsampled like
            % QUALITY_BITMASK and L2_QUALITY_BITMASK with the difference that
            % the final value could ADDITIONALLY POTENTIALLY later (in the
            % execution) be modified due to non-quality variable data
            % processing. The code is designed to anticipate this possibility.
            % Ex: >=1 timestamps/bin but below threshold.
            %       ==> mean=mstd=NaN.
            %       ==> QUALITY_FLAG=0? fill value?
            % Ex: >1 timestamps/bin and over threshold, some data samples=NaN
            %       ==> mean=mstd=NaN.
            %       ==> QUALITY_FLAG=0? fill value?
            %
            %
            % Correct zVar data types
            % -----------------------
            % "QUALITY_FLAG shall be a CDF_UINT1 flag"
            % "QUALITY_BITMASK shall be a CDF_UINT2 flag"
            % Source: SOL-SGS-TN-0009, "Metadata Definition for Solar Orbiter
            % Science Data"
            % --
            % "The optional CDF_UINT2 zVariable L2_QUALITY_BITMASK /.../"
            % Source:
            % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
            %====================================================================
            % Pre-allocate
            QUALITY_FLAG_dwns       = zeros(nRecordsDwns, 1, 'uint8');
            QUALITY_BITMASK_dwns    = zeros(nRecordsDwns, 1, 'uint16');
            L2_QUALITY_BITMASK_dwns = zeros(nRecordsDwns, 1, 'uint16');
            for i = 1:nRecordsDwns
                k = iRecordsInBinCa{i};

                QUALITY_FLAG_dwns(i) = ...
                    bicas.proc_sub23.downsample_bin_QUALITY_FLAG(...
                        InLfrCwfZv.QUALITY_FLAG( k) );

                % IMPLEMENTATION NOTE:
                % 2020-11-23: L2 zVar "QUALITY_BITMASK" is mistakenly
                % uint8/CDF_UINT1 when it should be uint16/CDF_UINT2. Must
                % therefore TYPECAST.
                % 2021-05-06: L2 zVar "QUALITY_BITMASK" type was changed from
                %   SOLO_L2_RPW-LFR-SURV-CWF-E_V11.skt: CDF_UTIN1
                % to 
                %   SOLO_L2_RPW-LFR-SURV-CWF-E_V12.skt: CDF_UINT2
                %   (SKELETON_MODS: V12=Feb 2021)
                % .
                QUALITY_BITMASK_dwns(i)    = ...
                    bicas.proc_sub23.downsample_bin_L12_QUALITY_BITMASK(...
                        uint16(...
                            InLfrCwfZv.QUALITY_BITMASK( k ))...
                        );

                L2_QUALITY_BITMASK_dwns(i) = ...
                    bicas.proc_sub23.downsample_bin_L12_QUALITY_BITMASK(...
                        InLfrCwfZv.L2_QUALITY_BITMASK( k ) );
            end



            %============================================================
            % Shared zVariables between all DOWNSAMPLED datasets
            %
            % (Initial value for QUALITY_FLAG; might be modified later.)
            %============================================================
            InitialDwnsZv = struct();
            InitialDwnsZv.Epoch              = zvEpochDwns;
            InitialDwnsZv.QUALITY_FLAG       = QUALITY_FLAG_dwns;
            InitialDwnsZv.QUALITY_BITMASK    = QUALITY_BITMASK_dwns;
            InitialDwnsZv.L2_QUALITY_BITMASK = L2_QUALITY_BITMASK_dwns;
            %
            % NOTE: Takes leap seconds into account.
            % NOTE/BUG: DELTA_PLUS_MINUS not perfect since the bin timestamp is
            % not centered for leap seconds. Epoch+-DELTA_PLUS_MINUS will thus
            % go outside/inside the bin boundaries for leap seconds. The same
            % problem exists for both positive and negative leap seconds.
            InitialDwnsZv.DELTA_PLUS_MINUS   = double(binSizeArrayNs / 2);
            
        end



        % Derive QUALITY_FLAG for ONE downsampled CDF record, from the
        % corresponding non-downsampled records (bin).
        %
        % NOTE: Handles empty bins.
        %
        function QUALITY_FLAG = downsample_bin_QUALITY_FLAG(zv_QUALITY_FLAG_segment)
            % TODO-DEC: Return NaN/fill value or 0 for empty bin?

            % IMPLEMENTATION NOTE: Just using min([zv_QUALITY_FLAG_segment; 0])
            % does not work.
            if isempty(zv_QUALITY_FLAG_segment)
                QUALITY_FLAG = 0;
            else
                QUALITY_FLAG = min(zv_QUALITY_FLAG_segment);
            end
        end



        % Derive a quality bitmask for ONE downsampled CDF record, from the
        % corresponding non-downsampled records (bin).
        %
        % NOTE: "L12_QUALITY_BITMASK" refers to both zVariables
        %   (1) QUALITY_BITMASK (set in L1), and
        %   (2) L2_QUALITY_BITMASK.
        %
        % NOTE: Handles empty bins.
        %
        function L12_QUALITY_BITMASK = downsample_bin_L12_QUALITY_BITMASK(...
                zv_L12_QUALITY_BITMASK_segment)
            % TODO-DEC: Return NaN/fill value or 0 for empty bin?

            % IMPLEMENTATION NOTE: 2020-11-23: L2 zVar "QUALITY_BITMASK" is
            % mistakenly uint8/CDF_UINT1 when it should be uint16/CDF_UINT2.
            assert(isa(zv_L12_QUALITY_BITMASK_segment, 'uint16'))

            if isempty(zv_L12_QUALITY_BITMASK_segment)
                L12_QUALITY_BITMASK = 0;
            else
                L12_QUALITY_BITMASK = bicas.utils.bitops.or(...
                    zv_L12_QUALITY_BITMASK_segment);
            end
        end

        

    end    % methods(Static, Access=private)

end
