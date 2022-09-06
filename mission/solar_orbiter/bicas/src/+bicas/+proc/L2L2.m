%
% Code for processing functions L2-->L2, i.e. downsampling
% SOLO_L2_RPW-LFR-SURV-CWF-E.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-18
%
classdef L2L2    
    % PROPOSAL: Automatic test code.
    
    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess = immutable)
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        
        
        
        % IMPLEMENTATION NOTE: It is not obvious whether this processing should
        % be run as a part of a separate s/w mode
        %   SOLO_L2_LFR-RPW-CWF-E --> DWNS,
        % or as part of the s/w mode
        %   SOLO_L1/L1R_LFR-RPW-CWF --> SOLO_L2_LFR-RPW-CWF-E (+DWNS).
        % The code is therefore designed so that it is easy to switch between
        % the two.
        %
        % IMPLEMENTATION NOTE: It is tempting to merge this function with
        % process_L2_to_L3() and have the same S/W MODE use it, since it
        % (1) has the same one input dataset, and
        % (2) produces a similarily downsampled dataset.
        % However, that might be a bad idea since it also
        % (1) uses another sampling rate (less shared processing), and
        % (2) should remain inofficial (both s/w mode and the output dataset),
        %     as opposed to process_L2_to_L3() which produces official datasets
        %     and might one day be "officially" run at ROC.
        %
        function OutLfrCwfDwns = process_LFRCWF_to_DWNS(InLfrCwf, SETTINGS, L)
            %
            % PROBLEM: How handle leap seconds if bin size <= 1 s?
            %   NOTE: Positive leap seconds are not a problem.
            %   PROPOSAL: Split bins WITH leap seconds? Then there is no
            %             problem(?).
            %
            % TODO-DEC: S/w mode?
            %   PROPOSAL: Separate mode
            %             SOLO_L2_LFR-RPW-CWF-E --> DWNS
            %       PRO: Can process files separately ==> faster.
            %       CON: Two separate modes require one SOLO_L2_LFR-RPW-CWF-E
            %            dataset.
            %           PRO: Can not automatically identify exactly one
            %               s/w mode from a set of input datasets as BPTD (external
            %               to BICAS) assumes for .
            %           PRO?!!: Has consequences for batch processing code?
            %
            %   PROPOSAL: Merged L1/L1R-->L2 mode 
            %             SOLO_L1/L1R_LFR-RPW-CWF --> SOLO_L2_LFR-RPW-CWF-E (+DWNS).
            %       PRO: Fewer modes for batch processing code to handle.
            %            ==> Faster identification of files+modes.
            %       PRO: More automatic synching of ORIS and DWNS versions.
            
            tTicToc = tic();



            %===========
            % Constants
            %===========
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            BIN_LENGTH_WOLS_NS        = int64(1e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);

            % 2021-05-24, YK: Only want to use QUALITY_FLAG>=2 data.
            QUALITY_FLAG_minForUse = SETTINGS.get_fv(...
                'PROCESSING.L2-CWF-DWNS.ZV_QUALITY_FLAG_MIN');


            
            %=========================
            % ~Generic initialization
            %=========================
            Ga = struct(...
                'OBS_ID',    InLfrCwf.Ga.OBS_ID, ...
                'SOOP_TYPE', InLfrCwf.Ga.SOOP_TYPE);
            %
            [Zv, iRecordsInBinCa] = bicas.proc.dwns.init_shared_downsampled(...
                InLfrCwf, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS, ...
                L);
            OutLfrCwfDwns = struct(...
                'Ga', Ga, ...
                'Zv', Zv);
            nRecordsDwns = numel(iRecordsInBinCa);
            
            
            
            zv_VDC = InLfrCwf.Zv.VDC;
            zv_EDC = InLfrCwf.Zv.EDC;
            nRecordsOris = numel(InLfrCwf.Zv.Epoch);
            
            
            
            % NOTE: Comparison will technically fail for QUALITY_FLAG fill
            % value, but that is acceptable (ideal result is ambiguous anyway).
            bDoNotUse = InLfrCwf.Zv.QUALITY_FLAG < QUALITY_FLAG_minForUse;
            
            zv_VDC(bDoNotUse, :) = NaN;
            zv_EDC(bDoNotUse, :) = NaN;
            
            
            
            clear InLfrCwf
            
            
            
            %===============================================
            % Downsample
            % ----------
            % NOTE: Exclude EAC, IBIAS1/2/3. /YK 2021-05-11
            %===============================================
            [OutLfrCwfDwns.Zv.VDC, ...
             OutLfrCwfDwns.Zv.VDCSTD] = bicas.proc.dwns.downsample_sci_zVar(...
                zv_VDC, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            
            [OutLfrCwfDwns.Zv.EDC, ...
             OutLfrCwfDwns.Zv.EDCSTD] = bicas.proc.dwns.downsample_sci_zVar(...
                zv_EDC, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            
            
            
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L2.process_LFRCWF_to_DWNS', tTicToc, ...
                nRecordsOris, 'ORIS record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L2.process_LFRCWF_to_DWNS', tTicToc, ...
                nRecordsDwns, 'DWNS record')

        end    % process_LFRCWF_to_DWNS



    end    % methods(Static)



end
