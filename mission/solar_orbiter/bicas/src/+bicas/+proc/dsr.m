%
% Collection of shared code used for creating downsampled datasets.
%
%
% NOTE
% ====
% 2020-11-23: L2 zVar "QUALITY_BITMASK" is
%   mistakenly uint8/CDF_UINT1 when it should be uint16/CDF_UINT2.
% 2023-08-30: Can neither see this for L1:QUALITY_BITMASK nor
%   L2:L2_QUALITY_BITMASK in
%   solo_L1_rpw-lfr-surv-cwf-cdag_20200212_V10.cdf
%   solo_L2_rpw-lfr-surv-cwf-e-cdag_20200229_V13.cdf
%   solo_L2_rpw-lfr-surv-cwf-e-cdag_20230823_V02.cdf
%   solo_L1_rpw-lfr-surv-cwf-cdag_20230823_V02.cdf
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-18
%
classdef dsr    
    %
    % PROPOSAL: Test code for code that downsamples.
    %   Ex: bicas.proc.dsr.downsample_ZV_bitmask()
    %   Ex: bicas.proc.dsr.downsample_ZV_minimum()
    %   --
    %   PRO: Can verify now uncertain edge cases.
    %           Ex: Quality ZVs when science data=fill values.
    %       PRO: Can verify bugfix for integer quality zVar=fill value when
    %            there is no science data.
    %
    % PROPOSAL: Work with FPA.
    %   PRO: Does not need to know fill value.
    %   CON: One FPA per bin? Slow?!
    %
    % PROPOSAL: Include bicas.utils.get_bin_indices() ?!
    %   CON: Potentially generic outside of BICAS.
    %       CON: Some code in this class can be considered similarly generic.
    %
    % PROPOSAL: Merge downsample_ZV_minimum() and downsample_ZV_bitmask()
    %           somwhow.
    %   PRO: Are very similar. Difference is the inner function called.



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Derive values which are used by all DOWNSAMPLED datasets.
        %
        %
        % RETURN VALUES
        % =============
        % InitialDsrZv
        %       Struct with zVariables. 
        % iRecordsInBinCa
        %       Distribution of non-downsampled records in bins.
        %
        function [InitialDsrZv, iRecordsInBinCa] = init_shared_DSR_ZVs(...
                InLfrCwfOsr, binLengthWolsNs, binTimestampPosWolsNs, L)
            
            % NOTE: Function argument InLfrCwfOsr contains too much information!
            % PROPOSAL: Only take argument for the needed variables.
            %   PROPOSAL: Wait until only using FPAs.
            %       PRO: Can abolish .ZvFv.
            
            tTicToc = tic();



            assert(isscalar(binLengthWolsNs))
            assert(isscalar(binTimestampPosWolsNs))
            
            %================
            % Calculate bins
            %================
            % Find bin boundary reference timestamp. This is used for
            % setting the bin boundaries together with the bin length.
            v = spdfbreakdowntt2000(InLfrCwfOsr.Zv.Epoch(1));
            % UTC subsecond (milliseconds, microseconds, nanoseconds)
            v(6)   = 5;   % UTC second
            v(7:9) = 0;   % Milliseconds, microseconds, nanoseconds
            boundaryRefTt2000 = spdfcomputett2000(v);
            % Find
            % (1) bin timestamps (downsampled timestamps to represent each bin),
            %     and
            % (2) which (non-downsampled) records belong to which bins
            %     (=downsampled records).
            [zvEpochDsr, iRecordsInBinCa, binSizeArrayNs] = ...
                bicas.proc.dsr.get_downsampling_bins(...
                    InLfrCwfOsr.Zv.Epoch, ...
                    boundaryRefTt2000, ...
                    binLengthWolsNs, ...
                    binTimestampPosWolsNs, ...
                    L);
            nRecordsOsr = numel(InLfrCwfOsr.Zv.Epoch);
            nRecordsDsr = numel(iRecordsInBinCa);
            
            
            
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
            % IMPLEMENTATION NOTE:
            % 2020-11-23: L2 zVar "L2_QUALITY_BITMASK" is mistakenly
            % uint8/CDF_UINT1 when it should be uint16/CDF_UINT2. Must
            % therefore TYPECAST.
            % 2021-05-06: L2 zVar "QUALITY_BITMASK" type was changed from
            %   SOLO_L2_RPW-LFR-SURV-CWF-E_V11.skt: CDF_UTIN1
            % to 
            %   SOLO_L2_RPW-LFR-SURV-CWF-E_V12.skt: CDF_UINT2
            %   (SKELETON_MODS: V12=Feb 2021)
            % .
            
%             tt = tic();
            zv_QUALITY_FLAG_FpaDsr    = bicas.proc.dsr.downsample_ZV_minimum(...
                InLfrCwfOsr.ZvFpa.QUALITY_FLAG,    iRecordsInBinCa);
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:QUALITY_FLAG', tt, ...
%                 nRecordsOsr, 'OSR record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:QUALITY_FLAG', tt, ...
%                 nRecordsDsr, 'DSR record')

%             tt = tic();
            zv_QUALITY_BITMASK_FpaDsr = bicas.proc.dsr.downsample_ZV_bitmask(...
                InLfrCwfOsr.ZvFpa.QUALITY_BITMASK, iRecordsInBinCa);
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:QUALITY_BITMASK', tt, ...
%                 nRecordsOsr, 'OSR record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:QUALITY_BITMASK', tt, ...
%                 nRecordsDsr, 'DSR record')

%             tt = tic();
            zv_L2_QUALITY_BITMASK_FpaDsr = bicas.proc.dsr.downsample_ZV_bitmask(...
                InLfrCwfOsr.ZvFpa.L2_QUALITY_BITMASK, iRecordsInBinCa);
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:L2_QUALITY_BITMASK', tt, ...
%                 nRecordsOsr, 'OSR record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.init_shared_DSR_ZVs:L2_QUALITY_BITMASK', tt, ...
%                 nRecordsDsr, 'DSR record')
            
            %============================================================
            % Shared zVariables between all DOWNSAMPLED datasets
            %
            % (Initial value for QUALITY_FLAG; might be modified later.)
            %============================================================
            InitialDsrZv = struct();
            InitialDsrZv.Epoch              = zvEpochDsr;
            InitialDsrZv.QUALITY_FLAG       = zv_QUALITY_FLAG_FpaDsr;
            InitialDsrZv.QUALITY_BITMASK    = zv_QUALITY_BITMASK_FpaDsr;
            InitialDsrZv.L2_QUALITY_BITMASK = zv_L2_QUALITY_BITMASK_FpaDsr;
            %
            % NOTE: Takes leap seconds into account.
            % NOTE/BUG: DELTA_PLUS_MINUS not perfect since the bin timestamp is
            % not centered for leap seconds. Epoch+-DELTA_PLUS_MINUS will thus
            % go outside/inside the bin boundaries for leap seconds. The same
            % problem exists for both positive and negative leap seconds.
            InitialDsrZv.DELTA_PLUS_MINUS   = double(binSizeArrayNs / 2);



            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.init_shared_DSR_ZVs', tTicToc, ...
                nRecordsOsr, 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.init_shared_DSR_ZVs', tTicToc, ...
                nRecordsDsr, 'DSR record')
        end
        
        
        
        % Utility function to help downsampling data by grouping together
        % adjacent time intervals that have the same length when discounting
        % leap seconds.
        %
        % NOTE: Function is only public so that automated test code can access
        % it.
        %
        %
        % ARGUMENTS
        % =========
        % zvAllTt2000
        %       Column array. ~Epoch.
        %       PROBLEM: Can not handle zvAllTt2000(1), zvAllTt2000(end)
        %       being during positive leap second.
        % boundaryRefTt2000
        %       Must not be during leap second. Specifies boundaries will be
        %       together with other arguments.
        % binLengthWolsNs
        %       Length of each bin.
        % binTimestampPosWolsNs
        %       Position of timestamp that represents bin, relative to beginning
        %       of bin.
        %
        %
        % RETURN VALUES
        % =============
        % zvTt2000
        %       Column array. Downsampled Epoch-like zVar. One timestamp per
        %       bin.
        % iRecordsInBinCa
        %       Column cell array.
        %       Indices to CDF records for respective bins.
        %       {iBin,1}(iSample,1) = Non-downsampled CDF record number.
        % binSizeArrayNs
        %       (iBin, 1). Bin size.
        %       RATIONALE: Useful for automatic testing, setting zVar
        %       DELTA_PLUS_MINUS (if one wants to account for leap seconds).
        %
        %
        % NAMING CONVENTIONS
        % ==================
        % See readme.txt
        % bin      : Time interval within which all corresponding CDF records
        %            should be condensed to one.
        % boundary : Edge of bin(s).
        %
        function [zvBinsTt2000, iRecordsInBinCa, binSizeArrayNs] = ...
            get_downsampling_bins(...
                zvAllTt2000, boundaryRefTt2000, ...
                binLengthWolsNs, binTimestampPosWolsNs, L)
            
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
            %
            % PROPOSAL: Better name that implies that it is not a generic
            %           function (works on Epoch with leap seconds).
            
            
            
            tTicToc = tic();
            

            
            % ASSERTIONS
            % ----------
            % Does not assert monotonic increase.
            bicas.utils.assert_ZV_Epoch(zvAllTt2000)
            % NOTE: Function algorithm assumes this monotonic increase.
            assert(issorted(zvAllTt2000, 'strictascend'))
            %
            bicas.utils.assert_ZV_Epoch(boundaryRefTt2000)
            assert(isscalar(boundaryRefTt2000))
            assert(isscalar(binLengthWolsNs))
            assert(isa(binLengthWolsNs,       'int64'))
            assert(isa(binTimestampPosWolsNs, 'int64'))
            assert((0 <= binTimestampPosWolsNs)...
                && (binTimestampPosWolsNs <= binLengthWolsNs))
            
            
            
            if isempty(zvAllTt2000)
                % CASE: zvAllTt2000 is empty.
                
                % NOTE: Later calls to
                % irf.cdf.time.TT2000_to_TT2000WOLS() are not applicable
                % if there is no timestamp. Must therefore have special case.
                
                zvBinsTt2000    = int64(ones(0,1));
                iRecordsInBinCa = cell(0,1);
                binSizeArrayNs  = int64(zeros(0,1));
                return
            end
            % CASE: zvAllTt2000 is not empty.
            
            
            
            ttw1           = irf.cdf.time.TT2000_to_TT2000WOLS(zvAllTt2000(1));
            ttw2           = irf.cdf.time.TT2000_to_TT2000WOLS(zvAllTt2000(end));
            boundaryRefTtw = irf.cdf.time.TT2000_to_TT2000WOLS(boundaryRefTt2000);
            
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
            
            boundariesTt2000 = irf.cdf.time.TT2000WOLS_to_TT2000(boundariesTtw);
            zvBinsTt2000     = irf.cdf.time.TT2000WOLS_to_TT2000(zvBinsTtw);
            binSizeArrayNs   = diff(boundariesTt2000);
            
            %===================================================================
            % Assign iRecordCa
            % ----------------
            % For every DSR record/bin, derive indices to the corresponding
            % OSR CDF records.
            %===================================================================
            % NOTE: DERIVING iRecordsInBinCa IS A SUBSTANTIAL PART OF THE
            % PROCESSING TIME (ESPECIALLY FOR LARGER DATASETS).
            iRecordsInBinCa = bicas.utils.get_bin_indices(...
                zvAllTt2000, boundariesTt2000, 10);



            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.get_downsampling_bins', tTicToc, ...
                numel(zvAllTt2000), 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.get_downsampling_bins', tTicToc, ...
                nBins,              'DSR record')
        end
        
        
        
        % Downsample a single NxM science zVar.
        %
        % Use bins and for every bin, derive median and MSTD over dimension 1
        % (within bin). Construct two zVariables for median+MSTD for the
        % corresponding downsampled CDF record.
        %
        % NOTE: Can handle zero records in bin.
        % NOTE: Function is only public so that automated test code can access
        % it.
        % NOTE: Can handle varying number of NaN in different zVar columns
        % (different channels).
        %
        % IMPLEMENTATION NOTE: Must count non-NaN samples for every channel (in
        % each bin) separately. This is a relevant use case since e.g. E-field
        % has Ex=NaN only, i.e. one can have bins with different number of NaN
        % in different columns.
        %
        %
        % ARGUMENTS
        % =========
        % zv
        %       (iCdfRecord, iChannel). Float.
        % nMinReqSamples
        %       Minimum number of samples (fill value or not) for not returning
        %       fill value.
        % iRecordsInBinCa
        %       Column cell array. {iBin, 1}
        %
        %
        % RETURN VALUES
        % =============
        % zvMed
        %       (iBin, iChannel). Median.
        % zvMstd
        %       (iBin, iChannel). MSTD.
        %
        function [zvMed, zvMstd] = downsample_sci_ZV(...
                zv, nMinReqRecords, iRecordsInBinCa, L)
            
            % PROPOSAL: Require nMinReqSamples >= 1? Code can handle 0, though it gives NaN.
            
            tTicToc = tic();
            
            

            % ASSERTION
            assert(isfloat(zv))
            assert(nMinReqRecords >= 0)
            [nRecordsOsr, nRecordsDsr, nSpr] = irf.assert.sizes(...
                zv,              [-1, -3], ...
                iRecordsInBinCa, [-2]);
            
            
            
            % Pre-allocate
            zvMed  = NaN(nRecordsDsr, nSpr);
            zvMstd = NaN(nRecordsDsr, nSpr);

            for iBin = 1:nRecordsDsr
                k           = iRecordsInBinCa{iBin};
                
                binZv       = zv(k, :);
                % Number of OSR records in bin.
                nBinRecords = size(binZv, 1);

                % (1) Pre-allocate (e.g. in case of zero samples/record).
                % (2) Set default values in case of too few records.
                binMed  = NaN(1, nSpr);
                binMstd = binMed;

                if nBinRecords >= nMinReqRecords
                    % CASE: Enough records, but not necessarily non-NaN samples.
                    
                    for iChannel = 1:nSpr
                        
                        % Vector of non-NaN samples, for one specific bin
                        % column.
                        tempSamples = binZv(:, iChannel);
                        samples     = tempSamples(~isnan(tempSamples));
                        
                        if size(samples, 1) < nMinReqRecords
                            % CASE: Too few non-NaN samples.
                            binMed(iChannel)  = NaN;
                            binMstd(iChannel) = NaN;
                        else
                            % CASE: Enough non-NaN samples.
                            binMed(iChannel)  = median(samples, 1);
                            binMstd(iChannel) = bicas.utils.mstd(...
                                samples, ...
                                binMed(iChannel), ...
                                1);
                        end
                    end
                end
                
                zvMed (iBin, :) = binMed;
                zvMstd(iBin, :) = binMstd;
                
                % clear binMed binMstd
            end    % for
            
            
            
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.downsample_sci_ZV', tTicToc, ...
%                 nRecordsOsr, 'OSR record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dsr.downsample_sci_ZV', tTicToc, ...
%                 nRecordsDsr,              'DSR record')
            
        end    % downsample_sci_ZV
        
        
        
        % ######################################################################
        % IMPLEMENTATIONS USING FPAs ALSO FOR BINS INTERNALLY
        % ######################################################################
        


%         % General-purpose function for downsampling data using bins.
%         % The implementation sends and receives FPAs to/from the inner function.
%         % Therefore tends to be very slow.
%         %
%         % ARGUMENTS
%         % =========
%         % ZvOsrFpa
%         %       FPA. ZV-like OSR data. 2D i.e. may have multiple columns.
%         % iRecordsInBinCa
%         % fhBin
%         %       Function handle. Condenses one bin of OSR CDF records to one (or
%         %       zero) CDF records.
%         %       ZvDsrFpa = fhBin(BinZvOsrFpa, FpFpa)
%         function DsrFpa = downsample_FPA(OsrFpa, iRecordsInBinCa, fhBin)
%             assert(isa(OsrFpa, 'bicas.utils.FillPositionsArray'), ...
%                 'Argument is not an instance of bicas.utils.FillPositionsArray.')
%             assert(ismatrix(OsrFpa))
%             assert(iscell(iRecordsInBinCa))            
%             assert(iscolumn(iRecordsInBinCa))
%             
%             nRecordsDsr = numel(iRecordsInBinCa);
%             
%             EmptyFpa = OsrFpa(1:0, :);    % 0x1;
%             FpFpa    = bicas.utils.FillPositionsArray.get_scalar_FP(OsrFpa.mc);
%             
%             dsrFpaCa = cell(nRecordsDsr, 1);   % Preallocate
% 
%             for iBin = 1:nRecordsDsr
%                 dsrFpaCa{iBin} = fhBin(...
%                     OsrFpa(iRecordsInBinCa{iBin}), FpFpa);
%             end
% 
%             DsrFpa = cat(1, EmptyFpa, dsrFpaCa{:});
%         end
%         
%         
%         
%         % Downsample a zVariable using pre-defined bins to minimum in each bin.
%         %
%         function DsrFpa = downsample_ZV_minimum_W_FPAs( OsrFpa, iRecordsInBinCa )
%             
%             function DsrFpa = bin_algo(BinOsrFpa, FpFpa)
%                 if isempty(BinOsrFpa)
%                     % CASE: Bin contains no values.
%                     DsrFpa = FpFpa;
%                 else
%                     % CASE: Bin contains non-zero number of values, which
%                     %       might be FPs.
% 
%                     nfpAr = BinOsrFpa.get_non_FP_data();
%                     if isempty(nfpAr)
%                         DsrFpa = FpFpa;
%                     else
%                         % NOTE: min(X, [], iDim)
%                         m        = min(nfpAr, [], 1);
%                         DsrFpa = bicas.utils.FillPositionsArray(m, 'NO_FILL_POSITIONS');
%                     end
%                 end
%             end
%             
%             assert(iscolumn(OsrFpa))
%             DsrFpa = bicas.proc.dsr.downsample_FPA(OsrFpa, iRecordsInBinCa, @bin_algo);
%         end
%         
%         
%         
%         % Downsample a zVariable using pre-defined bins to the logical OR of
%         % each bin.
%         %
%         function DsrFpa = downsample_ZV_bitmask_W_FPAs(OsrFpa, iRecordsInBinCa)
% 
%             function DsrFpa = bin_algo(BinOsrFpa, FpFpa)
%                 if isempty(BinOsrFpa)
%                     % CASE: Bin contains no values.
%                     DsrFpa = FpFpa;
%                 else
%                     % CASE: Bin contains non-zero number of values, which
%                     %       might be FPs.
% 
%                     nfpAr = BinOsrFpa.get_non_FP_data();
%                     if isempty(nfpAr)
%                         DsrFpa = FpFpa;
%                     else
%                         r = bicas.utils.bitops.or(nfpAr);
%                         DsrFpa = bicas.utils.FillPositionsArray(r, 'NO_FILL_POSITIONS');
%                     end
%                 end
%             end
% 
%             assert(iscolumn(OsrFpa))
%             DsrFpa = bicas.proc.dsr.downsample_FPA(OsrFpa, iRecordsInBinCa, @bin_algo);
%         end



        % ######################################################################
        % IMPLEMENTATIONS USING ARRAYS FOR BINS INTERNALLY
        % ######################################################################
        


        % General-purpose function for downsampling data using bins.
        %
        % ARGUMENTS
        % =========
        % OsrFpa
        %       FPA. ZV-like OSR data. 2D i.e. may have multiple columns.
        % iRecordsInBinCa
        % fhBin
        %       Function handle. Condenses one bin of OSR CDF records to one (or
        %       zero) CDF records.
        %       [binSamplesDsrAr, binFpAr] = fhBin(samplesOsrAr(iAr), fpOsrAr(iAr));
        function DsrFpa = downsample_W_INNER_ARRAYS(OsrFpa, iRecordsInBinCa, fhBin)
            assert(isa(OsrFpa, 'bicas.utils.FillPositionsArray'), ...
                'Argument is not an instance of bicas.utils.FillPositionsArray.')
            assert(ismatrix(OsrFpa))
            assert(iscell(iRecordsInBinCa))            
            assert(iscolumn(iRecordsInBinCa))
            
            nRecordsDsr = numel(iRecordsInBinCa);
            
            samplesOsrAr   = OsrFpa.get_data();
            fpOsrAr        = OsrFpa.fpAr;
            
            samplesDsrArCa = cell(nRecordsDsr, 1);    % Preallocate
            fpDsrArCa      = cell(nRecordsDsr, 1);    % Preallocate

            for iBin = 1:nRecordsDsr
                iAr = iRecordsInBinCa{iBin};

                [binSamplesDsrAr, binFpAr] = fhBin(samplesOsrAr(iAr), fpOsrAr(iAr));
                
                if 1
                    % DEBUG: Verify the return values from function.
                    assert(isa(      binSamplesDsrAr, OsrFpa.mc))
                    assert(islogical(binFpAr        ))
                    assert(isequaln(...
                        size(binSamplesDsrAr), ...
                        size(binFpAr        )))
                end
                
                samplesDsrArCa{iBin} = binSamplesDsrAr;
                fpDsrArCa{     iBin} = binFpAr;
            end

            samplesDsrAr = cat(1, zeros(0, 1, OsrFpa.mc), samplesDsrArCa{:});
            fpDsrAr      = cat(1, false(0, 1),            fpDsrArCa{     :});
            
            DsrFpa = bicas.utils.FillPositionsArray(samplesDsrAr, 'FILL_POSITIONS', fpDsrAr);
        end



        % Downsample a zVariable using pre-defined bins to minimum in each bin.
        %
        function [DsrFpa] = downsample_ZV_minimum( OsrFpa, iRecordsInBinCa )
           
            function [binSamplesDsrAr, binFpDsrAr] = bin_algo(binSamplesOsrAr, binFpOsrAr)

                binSamplesOsrAr(binFpOsrAr, :) = [];   % Delete rows with FPs.

                if isempty(binSamplesOsrAr)
                    % CASE: Bin contains no values.
                    binSamplesDsrAr = fv;
                    binFpDsrAr      = true;
                else
                    % CASE: Bin contains non-zero number of values, which
                    %       might be FPs.

                    if isempty(binSamplesOsrAr)
                        binSamplesDsrAr = fv;
                        binFpDsrAr      = true;
                    else
                        % NOTE: min(X, [], iDim)
                        binSamplesDsrAr = min(binSamplesOsrAr, [], 1);
                        binFpDsrAr      = false;
                    end
                end
            end
            
            assert(iscolumn(OsrFpa))
            fv = zeros(1,1, OsrFpa.mc);   % Used by inner function.
            
            DsrFpa = bicas.proc.dsr.downsample_W_INNER_ARRAYS(OsrFpa, iRecordsInBinCa, @bin_algo);
        end


        
        % Downsample a zVariable using pre-defined bins to the logical OR of
        % each bin.
        %
        function DsrFpa = downsample_ZV_bitmask(OsrFpa, iRecordsInBinCa)

            function [binSamplesDsrAr, binFpDsrAr] = bin_algo(binSamplesOsrAr, binFpOsrAr)

                binSamplesOsrAr(binFpOsrAr, :) = [];   % Delete rows with FPs.

                if isempty(binSamplesOsrAr)
                    % CASE: Bin contains no values.
                    binSamplesDsrAr = fv;
                    binFpDsrAr      = true;
                else
                    % CASE: Bin contains non-zero number of values, which
                    %       might be FPs.

                    if isempty(binSamplesOsrAr)
                        binSamplesDsrAr = fv;
                        binFpDsrAr      = true;
                    else
                        binSamplesDsrAr = bicas.utils.bitops.or(binSamplesOsrAr);
                        binFpDsrAr      = false;
                    end
                end
            end

            assert(iscolumn(OsrFpa))
            fv = zeros(1,1, OsrFpa.mc);   % Used by inner function.

            DsrFpa = bicas.proc.dsr.downsample_W_INNER_ARRAYS(OsrFpa, iRecordsInBinCa, @bin_algo);
        end
        
        
        
    end


    
end
