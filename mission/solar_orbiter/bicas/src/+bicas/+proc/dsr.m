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



        % Derive template ZVs (not GA) for all datasets based on DOWNSAMPLED LFR
        % CWF.
        %
        %
        % RETURN VALUES
        % =============
        % TemplateDsrZv
        %       Struct with common ZVs:
        %           Epoch
        %           QUALITY_FLAG
        %           QUALITY_BITMASK
        %           L2_QUALITY_BITMASK
        %           DELTA_PLUS_MINUS
        % iRecordsInBinCa
        %       Distribution of non-downsampled records in bins.
        %
        function [TemplateDsrZv, iRecordsInBinCa] = get_LFR_CWF_DSR_ZVs_template(...
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
            
            zv_QUALITY_FLAG_FpaDsr    = bicas.proc.dsr.downsample_ZV_minimum(...
                InLfrCwfOsr.ZvFpa.QUALITY_FLAG,    iRecordsInBinCa);

            zv_QUALITY_BITMASK_FpaDsr = bicas.proc.dsr.downsample_ZV_bitmask(...
                InLfrCwfOsr.ZvFpa.QUALITY_BITMASK, iRecordsInBinCa);

            zv_L2_QUALITY_BITMASK_FpaDsr = bicas.proc.dsr.downsample_ZV_bitmask(...
                InLfrCwfOsr.ZvFpa.L2_QUALITY_BITMASK, iRecordsInBinCa);
            
            %============================================================
            % Shared zVariables between all DOWNSAMPLED datasets
            %
            % (Initial value for QUALITY_FLAG; might be modified later.)
            %============================================================
            Zv = struct();
            Zv.Epoch              = zvEpochDsr;
            Zv.QUALITY_FLAG       = zv_QUALITY_FLAG_FpaDsr;
            Zv.QUALITY_BITMASK    = zv_QUALITY_BITMASK_FpaDsr;
            Zv.L2_QUALITY_BITMASK = zv_L2_QUALITY_BITMASK_FpaDsr;
            %
            % NOTE: Takes leap seconds into account.
            % NOTE/BUG: DELTA_PLUS_MINUS not perfect since the bin timestamp is
            % not centered for leap seconds. Epoch+-DELTA_PLUS_MINUS will thus
            % go outside/inside the bin boundaries for leap seconds. The same
            % problem exists for both positive and negative leap seconds.
            Zv.DELTA_PLUS_MINUS   = bicas.utils.FPArray(binSizeArrayNs / 2);
            
            TemplateDsrZv = Zv;



            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template', tTicToc, ...
                nRecordsOsr, 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template', tTicToc, ...
                nRecordsDsr, 'DSR record')
        end
        
        
        
        % Utility function to help downsampling data by grouping together
        % adjacent time intervals that have the same length when discounting
        % leap seconds.
        %
        %
        % ARGUMENTS
        % =========
        % zvAllTt2000
        %       Column array. ~Epoch. Timestamps of OSR samples.
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
            %
            % PROPOSAL: Same-sized bins when including leap seconds.
            %   PRO: More "clean".
            %       PRO: Same-length bins.
            %       PRO: No problems with leap seconds(?).
            %       CON: Bin boundaries will not line up if merging data from
            %            multiple days.
            %
            % PROPOSAL: Separate function for generating boundaries.
            
            tTicToc = tic();
            

            
            % ASSERTIONS
            % ----------
            % Function does not assert strict increase.
            bicas.utils.assert_ZV_Epoch(zvAllTt2000)
            % NOTE: Function algorithm assumes this strict increase.
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
            % "Round" ttw1 down to nearest lower interval boundary, assuming
            % that boundaries are always on the form TTW=n*binLengthWolsNs.
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
        
        
        
        % Downsample a single NxM science ZV.
        %
        % Use bins and for every bin, derive median and MSTD over dimension 1
        % (within bin). Construct two zVariables for median+MSTD for the
        % corresponding downsampled CDF record.
        %
        % NOTE: Can handle zero records in bin.
        % NOTE: Function is only public so that automated test code can access
        %       it.
        % NOTE: Can handle varying number of FP in different zVar columns
        %       (different channels).
        %
        % IMPLEMENTATION NOTE: Must count NFP samples for every channel (in each
        % bin) separately. This is a relevant use case since e.g. E-field has
        % E_x=FV/FP only, i.e. one can have bins with different number of FP in
        % different columns.
        %
        % IMPLEMENTATION NOTE: Could generalize function to arbitrary number
        % of dimensions, but that is not needed. Therefore only requires 2D
        % data (iRecord, iChannel).
        %
        %
        % ARGUMENTS
        % =========
        % OsrFpa
        %       (iCdfRecord, iChannel). FPA double.
        % nMinNfpSamplesPerBin
        %       Minimum number of NFP samples (per channel/column) per bin for
        %       not returning FP.
        % iRecordsInBinCa
        %       Column cell array. {iBin, 1}
        %
        %
        % RETURN VALUES
        % =============
        % MedianDsrFpa
        %       (iBin, iChannel). FPA. Median.
        % MstdDsrFpa
        %       (iBin, iChannel). FPA. MSTD.
        %
        function [MedianDsrFpa, MstdDsrFpa] = downsample_sci_ZV(...
                OsrFpa, nMinNfpSamplesPerBin, iRecordsInBinCa, L)
            
            % PROPOSAL: Require nMinNfpSamplesPerBin >= 1? Code can handle 0, though it gives NaN.
            
            % IMPLEMENTATION NOTE: Not using bicas.proc.dsr.downsample() since
            % it is not designed for returning two values per bin.
            
            tTicToc = tic();

            % ASSERTIONS
            assert(isa(OsrFpa, 'bicas.utils.FPArray'))
            assert(strcmp(OsrFpa.mc, 'double'))
            assert(nMinNfpSamplesPerBin >= 0)
            [nRecordsOsr, nBins, nSpr] = irf.assert.sizes(...
                OsrFpa,          [-1, -3], ...
                iRecordsInBinCa, [-2]);
            nRecordsDsr = numel(iRecordsInBinCa);
            
            

            % FPA --> OSR arrays
            % ------------------
            % NOTE: Convert to double-NaN since that is what the algorithm
            % works with (MSTD, really).
            dataOsr = OsrFpa.array(NaN);
            fpOsr   = OsrFpa.fpAr();
            
            % Pre-allocate DSR arrays.
            medianDsr = zeros(nBins, nSpr, OsrFpa.mc);
            mstdDsr   = zeros(nBins, nSpr, OsrFpa.mc);
            
            % Pre-allocate DSR bin arrays.
            binDsrSize    = size(dataOsr);
            binDsrSize(1) = 1;
            binMedian = zeros(binDsrSize);
            binMstd   = zeros(binDsrSize);
            
            for iBin = 1:nBins
                % OSR bin arrays
                kBin = iRecordsInBinCa{iBin};
                binDataOsr = dataOsr(kBin, :);
                binFpOsr   = fpOsr  (kBin, :);

                for iChannel = 1:nSpr

                    % Vector of NFP samples (for one channel/bin column).
                    binChannelData    = binDataOsr              (:, iChannel);
                    binChannelDataNfp = binChannelData(~binFpOsr(:, iChannel));

                    if numel(binChannelDataNfp) < nMinNfpSamplesPerBin
                        % CASE: Too few NFP samples.
                        binMedian(1, iChannel) = NaN;
                        binMstd  (1, iChannel) = NaN;
                    else
                        % CASE: Enough NFP samples.
                        medianScalarData = median(binChannelDataNfp, 1);
                        mstdScalarData   = bicas.utils.mstd(...
                            binChannelDataNfp, ...
                            medianScalarData, ...
                            1);
                        
                        binMedian(1, iChannel) = medianScalarData;   % May be NaN.
                        binMstd  (1, iChannel) = mstdScalarData;     % May be NaN.
                    end
                end
                
                medianDsr(iBin, :) = binMedian;
                mstdDsr  (iBin, :) = binMstd;
            end    % for

            % NOTE: Must be (MATLAB class) double FPAs.
            MedianDsrFpa = bicas.utils.FPArray(medianDsr, 'FILL_VALUE', NaN);
            MstdDsrFpa   = bicas.utils.FPArray(mstdDsr,   'FILL_VALUE', NaN);



            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.downsample_sci_ZV', tTicToc, ...
                nRecordsOsr, 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.dsr.downsample_sci_ZV', tTicToc, ...
                nRecordsDsr,              'DSR record')
            
        end    % downsample_sci_ZV



        % General-purpose function for downsampling data using bins.
        %
        % ARGUMENTS
        % =========
        % OsrFpa
        %       FPA. ZV-like OSR data. 2D i.e. may have multiple columns.
        % iRecordsInBinCa
        % fhBinToRecord
        %       Function handle. Condenses one bin of OSR CDF records to exactly
        %       one CDF record.
        %       [binSamplesDsrAr, binFpAr] = fhBinToRecord(samplesOsrAr, fpOsrAr);
        function DsrFpa = downsample(OsrFpa, iRecordsInBinCa, fhBinToRecord)
            % PROPOSAL: Abolish
            
            assert(isa(OsrFpa, 'bicas.utils.FPArray'), ...
                'Argument is not an instance of bicas.utils.FPArray.')
            assert(ismatrix(OsrFpa))
            assert(iscell(iRecordsInBinCa))            
            assert(iscolumn(iRecordsInBinCa))
            
            nRecordsDsr = numel(iRecordsInBinCa);
            
            samplesOsrAr   = OsrFpa.array();
            fpOsrAr        = OsrFpa.fpAr;
            
            samplesDsrArCa = cell(nRecordsDsr, 1);    % Preallocate
            fpDsrArCa      = cell(nRecordsDsr, 1);    % Preallocate

            for iBin = 1:nRecordsDsr
                iAr = iRecordsInBinCa{iBin};

                [binSamplesDsrAr, binFpAr] = fhBinToRecord(samplesOsrAr(iAr), fpOsrAr(iAr));
                
%                 if 0
%                     % DEBUG: Verify the return values from function.
%                     assert(isa(      binSamplesDsrAr, OsrFpa.mc))
%                     assert(islogical(binFpAr        ))
%                     assert(size(binSamplesDsrAr, 1)  == 1)
%                     assert(isequaln(...
%                         size(binSamplesDsrAr), ...
%                         size(binFpAr        )))
%                 end
                
                samplesDsrArCa{iBin} = binSamplesDsrAr;
                fpDsrArCa{     iBin} = binFpAr;
            end

            samplesDsrAr = cat(1, zeros(0, 1, OsrFpa.mc), samplesDsrArCa{:});
            fpDsrAr      = cat(1, false(0, 1),            fpDsrArCa{     :});
            
            DsrFpa = bicas.utils.FPArray(samplesDsrAr, 'FILL_POSITIONS', fpDsrAr);
        end



        % Downsample a zVariable using pre-defined bins to minimum in each bin.
        %
        function [DsrFpa] = downsample_ZV_minimum( OsrFpa, iRecordsInBinCa )
            assert(isa(OsrFpa, 'bicas.utils.FPArray'))
           
            function [binSamplesDsrAr, binFpDsrAr] = bin_to_record(binSamplesOsrAr, binFpOsrAr)

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
            
            DsrFpa = bicas.proc.dsr.downsample(OsrFpa, iRecordsInBinCa, @bin_to_record);
        end


        
        % Downsample a zVariable (FPA) using pre-defined bins to the "logical OR
        % value" each OSR bin.
        %
        function DsrFpa = downsample_ZV_bitmask(OsrFpa, iRecordsInBinCa)
            assert(isa(OsrFpa, 'bicas.utils.FPArray'))

            % How to reduce one bin into one record.
            function [binSamplesDsrAr, binFpDsrAr] = bin_to_record(binSamplesOsrAr, binFpOsrAr)

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

            DsrFpa = bicas.proc.dsr.downsample(OsrFpa, iRecordsInBinCa, @bin_to_record);
        end
        
        
        
    end


    
end
