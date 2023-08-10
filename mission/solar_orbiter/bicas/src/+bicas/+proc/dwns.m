%
% Collection of shared code used for creating downsampled datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-18
%
classdef dwns    
    %
    % PROPOSAL: Test code for code that downsamples.
    %   Ex: bicas.proc.dwns.downsample_bin_L12_QUALITY_BITMASK() -- Too trivial?
    %   Ex: bicas.proc.dwns.downsample_bin_QUALITY_FLAG()        -- Too trivial?
    %   --
    %   PRO: Can verify now uncertain edge cases.
    %           Ex: Quality zVars when science data=fill values.
    %       PRO: Can verify bugfix for integer quality zVar=fill value when
    %            there is no science data.
    %
    % PROPOSAL: Include bicas.utils.get_bin_indices ?!
    %   CON: Potentially generic outside of BICAS.



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
        % InitialDwnsZv
        %       Struct with zVariables. 
        % iRecordsInBinCa
        %       Distribution of non-downsampled records in bins.
        %
        function [InitialDwnsZv, iRecordsInBinCa] = init_shared_downsampled(...
                InLfrCwf, binLengthWolsNs, binTimestampPosWolsNs, L)
            
            tTicToc = tic();



            assert(isscalar(binLengthWolsNs))
            assert(isscalar(binTimestampPosWolsNs))
            
            %================
            % Calculate bins
            %================
            % Find bin boundary reference timestamp. This is used for
            % setting the bin boundaries together with the bin length.
            v = spdfbreakdowntt2000(InLfrCwf.Zv.Epoch(1));
            % UTC subsecond (milliseconds, microseconds, nanoseconds)
            v(6)   = 5;   % UTC second
            v(7:9) = 0;
            boundaryRefTt2000 = spdfcomputett2000(v);
            % Find
            % (1) bin timestamps (downsampled timestamps to represent each bin),
            %     and
            % (2) which (non-downsampled) records belong to which bins
            %     (=downsampled records).
            [zvEpochDwns, iRecordsInBinCa, binSizeArrayNs] = ...
                bicas.proc.dwns.get_downsampling_bins(...
                    InLfrCwf.Zv.Epoch, ...
                    boundaryRefTt2000, ...
                    binLengthWolsNs, ...
                    binTimestampPosWolsNs, ...
                    L);
            nRecordsOris = numel(InLfrCwf.Zv.Epoch);
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
            for iBin = 1:nRecordsDwns
                k = iRecordsInBinCa{iBin};

                QUALITY_FLAG_dwns(iBin) = ...
                    bicas.proc.dwns.downsample_bin_QUALITY_FLAG(...
                        InLfrCwf.Zv.QUALITY_FLAG( k ), ...
                        InLfrCwf.ZvFv.QUALITY_FLAG);

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
                QUALITY_BITMASK_dwns(iBin)    = ...
                    bicas.proc.dwns.downsample_bin_L12_QUALITY_BITMASK(...
                        uint16( InLfrCwf.Zv.QUALITY_BITMASK( k ) ), ...
                        InLfrCwf.ZvFv.QUALITY_BITMASK);

                L2_QUALITY_BITMASK_dwns(iBin) = ...
                    bicas.proc.dwns.downsample_bin_L12_QUALITY_BITMASK(...
                        InLfrCwf.Zv.L2_QUALITY_BITMASK( k ), ...
                        InLfrCwf.ZvFv.L2_QUALITY_BITMASK);
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



%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dwns.init_shared_downsampled', tTicToc, ...
%                 nRecordsOris, 'ORIS record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dwns.init_shared_downsampled', tTicToc, ...
%                 nRecordsDwns, 'DWNS record')
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
        %       Column array. Downsampled Epoch-like zVar. One timestamp per
        %       bin.
        % iRecordsInBinCa
        %       Column cell array.
        %       Indices to CDF records for respective bins.
        %       {iBin, 1}(iSamples,1) = Non-downsampled CDF record number.
        % binSizeArrayNs
        %       (iBin, 1). Bin size.
        %       RATIONALE: Useful for automatic testing, setting zVar
        %       DELTA_PLUS_MINUS (if one wants to account for leap seconds).
        %
        %
        % NAMING CONVENTIONS
        % ==================
        % See readme.txt
        % bin      = Time interval within which all corresponding CDF records
        %            should be condensed to one.
        % boundary = Edge of bin(s).
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
            bicas.utils.assert_zv_Epoch(zvAllTt2000)
            % NOTE: Function algorithm assumes this monotonic increase.
            assert(issorted(zvAllTt2000, 'strictascend'))
            %
            bicas.utils.assert_zv_Epoch(boundaryRefTt2000)
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
            % For every DWNS record/bin, derive indices to the corresponding
            % ORIS CDF records.
            %===================================================================
            % NOTE: DERIVING iRecordsInBinCa IS A SUBSTANTIAL PART OF THE
            % PROCESSING TIME (ESPECIALLY FOR LARGER DATASETS).
            iRecordsInBinCa = bicas.utils.get_bin_indices(...
                zvAllTt2000, boundariesTt2000, 10);



            bicas.log_speed_profiling(L, ...
                'bicas.proc.dwns.get_downsampling_bins', tTicToc, ...
                numel(zvAllTt2000), 'ORIS record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.dwns.get_downsampling_bins', tTicToc, ...
                nBins,              'DWNS record')
        end
        
        
        
        % Downsample a single NxM science zVar.
        %
        % Use bins and for every bin, derive median and modified standard
        % deviation over dimension 1 (within bin). Construct two zVariables for
        % median+MSTD for the corresponding downsampled CDF record.
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
        % zvMed  : (iBin, iChannel). Median.
        % zvMstd : (iBin, iChannel). Modified STandard Deviation (MSTD).
        %
        function [zvMed, zvMstd] = downsample_sci_zVar(...
                zv, nMinReqRecords, iRecordsInBinCa, L)
            
            % PROPOSAL: Require nMinReqSamples >= 1? Code can handle 0, though it gives NaN.
            
            
            
            tTicToc = tic();
            
            

            % ASSERTION
            assert(isfloat(zv))
            assert(nMinReqRecords >= 0)
            [nRecordsOris, nRecordsDwns, nSpr] = irf.assert.sizes(...
                zv,              [-1, -3], ...
                iRecordsInBinCa, [-2]);
            
            
            
            % Pre-allocate
            zvMed  = NaN(nRecordsDwns, nSpr);
            zvMstd = NaN(nRecordsDwns, nSpr);

            for iBin = 1:nRecordsDwns
                k           = iRecordsInBinCa{iBin};
                
                binZv       = zv(k, :);
                % Number of ORIS records in bin.
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
                            binMstd(iChannel) = bicas.utils.modif_std_deviation(...
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
%                 'bicas.proc.dwns.downsample_sci_zVar', tTicToc, ...
%                 nRecordsOris, 'ORIS record')
%             bicas.log_speed_profiling(L, ...
%                 'bicas.proc.dwns.downsample_sci_zVar', tTicToc, ...
%                 nRecordsDwns,              'DWNS record')
            
        end    % downsample_sci_zVar



    end
        

    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Derive QUALITY_FLAG for ONE downsampled CDF record, from the
        % corresponding non-downsampled records (bin).
        %
        % NOTE: Handles empty bins.
        %
        function QUALITY_FLAG = downsample_bin_QUALITY_FLAG(...
                zv_QUALITY_FLAG_bin, fillValue)

            % Remove records with fill values.
            bUse = (zv_QUALITY_FLAG_bin ~= fillValue);
            zv_QUALITY_FLAG_bin = zv_QUALITY_FLAG_bin(bUse, :, :);
            
            % IMPLEMENTATION NOTE: Using min([zv_QUALITY_FLAG_segment; 0])
            % does not work if one wants to return 0 for empty bins.
            if isempty(zv_QUALITY_FLAG_bin)
                QUALITY_FLAG = fillValue;
            else
                % CASE: zv_QUALITY_FLAG_bin contains no fill values.
                QUALITY_FLAG = min(zv_QUALITY_FLAG_bin);
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
                zv_L12_QUALITY_BITMASK_bin, fillValue)

            % IMPLEMENTATION NOTE: 2020-11-23: L2 zVar "QUALITY_BITMASK" is
            % mistakenly uint8/CDF_UINT1 when it should be uint16/CDF_UINT2.
            assert(isa(zv_L12_QUALITY_BITMASK_bin, 'uint16'))
            
            % Remove records with fill values.
            b = zv_L12_QUALITY_BITMASK_bin ~= fillValue;
            zv_L12_QUALITY_BITMASK_bin = zv_L12_QUALITY_BITMASK_bin(b, :, :);

            if isempty(zv_L12_QUALITY_BITMASK_bin)
                L12_QUALITY_BITMASK = fillValue;
            else
                L12_QUALITY_BITMASK = bicas.utils.bitops.or(...
                    zv_L12_QUALITY_BITMASK_bin);
            end
        end

        

    end

    
    
end
