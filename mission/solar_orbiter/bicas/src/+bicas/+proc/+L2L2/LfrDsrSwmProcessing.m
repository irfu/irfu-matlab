%
% SWMP for downsampling SOLO_L2_RPW-LFR-SURV-CWF-E L2-->L2 (unofficial
% output datasets).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef LfrDsrSwmProcessing < bicas.proc.SwmProcessing
    % PROPOSAL: Automatic test code.



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)
        
        
        
        % OVERRIDE
        function OutputDatasetsMap = production_function(obj, ...
            InputDatasetsMap, rctDir, NsoTable, Bso, L)
            
            InLfrCwf = InputDatasetsMap('OSR_cdf');
            
            OutLfrCwfDsr = bicas.proc.L2L2.LfrDsrSwmProcessing.process_LFRCWF_to_DSR(InLfrCwf, Bso, L);
            
            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('DSR_cdf') = OutLfrCwfDsr;
        end



    end    % methods(Access=public)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % IMPLEMENTATION NOTE: It is not obvious whether this processing should
        % be run as a part of a separate s/w mode
        %   SOLO_L2_LFR-RPW-CWF-E --> DSR,
        % or as part of the s/w mode
        %   SOLO_L1/L1R_LFR-RPW-CWF --> SOLO_L2_LFR-RPW-CWF-E (+DSR).
        % The code is therefore designed so that it is easy to switch between
        % the two.
        %
        % IMPLEMENTATION NOTE: It is tempting to merge this function with
        % process_L2_to_L3() and have the same S/W MODE use it, since it
        % (1) has the same one input dataset, and
        % (2) produces a similarily downsampled dataset.
        % However, that might be a bad idea since it also
        % (1) uses another sampling rate (less shared processing), and
        % (2) should remain unofficial (both s/w mode and the output dataset),
        %     as opposed to process_L2_to_L3() which produces official datasets
        %     and might one day be "officially" run at ROC.
        %
        function OutLfrCwfDsr = process_LFRCWF_to_DSR(InLfrCwf, Bso, L)
            %
            % PROBLEM: How handle leap seconds if bin size <= 1 s?
            %   NOTE: Positive leap seconds are not a problem.
            %   PROPOSAL: Split bins WITH leap seconds? Then there is no
            %             problem(?).
            
            tTicToc = tic();



            %===========
            % Constants
            %===========
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            BIN_LENGTH_WOLS_NS        = int64(1e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);

            % 2021-05-24, YK: Only want to use QUALITY_FLAG>=2 data.
            QUALITY_FLAG_minForUse = uint8(Bso.get_fv(...
                'PROCESSING.L2-CWF-DSR.ZV_QUALITY_FLAG_MIN'));


            
            %=========================
            % ~Generic initialization
            %=========================
            Ga = struct(...
                'OBS_ID',    InLfrCwf.Ga.OBS_ID, ...
                'SOOP_TYPE', InLfrCwf.Ga.SOOP_TYPE);
            %
            [Zv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
                InLfrCwf, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS, ...
                L);
            OutLfrCwfDsr = struct(...
                'Ga', Ga, ...
                'Zv', Zv);
            nRecordsDsr = numel(iRecordsInBinCa);
            
            
            
            VdcOsrFpa = InLfrCwf.ZvFpa.VDC;
            EdcOsrFpa = InLfrCwf.ZvFpa.EDC;
            nRecordsOsr = numel(InLfrCwf.Zv.Epoch);
            
            
            
            % NOTE: Unclear how treat QUALITY_FLAG=FV.
            bDoNotUseFpa = InLfrCwf.ZvFpa.QUALITY_FLAG < QUALITY_FLAG_minForUse;
            bDoNotUse    = bDoNotUseFpa.array(false);   % FV = false wise?
            
            VdcOsrFpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;
            EdcOsrFpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;
            
            
            
            clear InLfrCwf
            
            
            
            %===============================================
            % Downsample
            % ----------
            % NOTE: Exclude EAC, IBIAS1/2/3. /YK 2021-05-11
            %===============================================
            [VdcDsrFpa, VdcstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                VdcOsrFpa.cast('double'), ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutLfrCwfDsr.Zv.VDC    = VdcDsrFpa.cast('single');
            OutLfrCwfDsr.Zv.VDCSTD = VdcstdDsrFpa.cast('single');
            
            [EdcDsrFpa, EdcstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                EdcOsrFpa.cast('double'), ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutLfrCwfDsr.Zv.EDC    = EdcDsrFpa.cast('single');
            OutLfrCwfDsr.Zv.EDCSTD = EdcstdDsrFpa.cast('single');

            

            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L2.process_LFRCWF_to_DSR', tTicToc, ...
                nRecordsOsr, 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L2.process_LFRCWF_to_DSR', tTicToc, ...
                nRecordsDsr, 'DSR record')

        end    % process_LFRCWF_to_DSR



    end    % methods(Static, Access=private)



end
