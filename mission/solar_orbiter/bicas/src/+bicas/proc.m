%
% Set of functions for producing one specific output dataset PDV from the necessary input dataset PDVs.
%
%
% PRODUCTION FUNCTIONS
% ====================
% A function with interface
%   OutputsMap = produce_*(InputsMap, Cal)
% where
% Cal        : A bicas.calib object.
% InputsMap  : containers.Map with
%                <keys>       : String defining a name of an input ("prodFuncInputKey" in swmode_defs).
%                <values>     : A struct with data corresponding to a CDF file (zVariables+global attributes).
% OutputsMap : containers.Map with
%                <keys>       : String defining a name of an output ("prodFuncOutputKey" in swmode_defs).
%                <values>     : A struct with data corresponding to a CDF file (zVariables).
% NOTE: In practice, anonymous functions with the correct interface are used to wrap the actual implementing functions
% (with another interface).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-07-30
%
classdef proc
    % PROPOSAL: Other name of class.
    %   PROPOSAL: production_functions
    %   PROPOSAL: production
    %   PROPOSAL: processing
    %       TODO-DECISION: Name relationship to proc_sub, proc_utils?
    %           PROPOSAL: Rename proc_sub
    %               processing_functions
    %               processing_subfunctions
    %           PROPOSAL: Rename proc_utils
    %               processing_utils
    % --
    % TODO-DECISION: Term and naming scheme for these "complete processing functions"?
    %   NOTE: Already used term "processing function" and process_* for partal processing.
    %   PROPOSAL: produce_*
    %   PROPOSAL: generate_*
    %   PROPOSAL: derive_*
    %   PROPOSAL: mode_*
    %   PROPOSAL: swmode_*
    %       CON: Already considering changing the term "s/w mode".
    %   PROPOSAL: pipeline_*
    %       CON: Conflicts with use of RODP, ROC-SGSE pipelines.
    %   PROPOSAL: swmode_pipelines
    %   PROPOSAL: process_*
    %       CON: Already used in proc_sub.
    %   TODO-DECISION: Include skeleton version in function names?
    % --
    % PROPOSAL: Production functions should not assume/specify any particular input dataset version, but read it out
    %           from global attributes (part of the PDV).
    % PROPOSAL: Somehow associate metadata with each function.
    %   PRO: For a given OUTPUT dataset/PDID, one can get the list of possible input datasets/PDIDs
    %   PROPOSAL: A given production function with a particular output dataset/PDID represents a s/w mode and will (for
    %       the right arguments, e.g. empty) return information on the corresponding inputs (incl. for s/w descriptor).
    % 
    %   NOTE: Metadata associated with each s/w mode

    
    
    methods(Static, Access=public)
        
        
        
        % ARGUMENTS
        % =========
        % InputDatasetsMap : containers.Map: key=<argument key> --> value=PDV
        %                    for input CDF.
        % inputSciDsi      : The science input dataset will be interpreted as
        %                    having this DATASET_ID.
        %                    RATIONALE: InputDatasetsMap should contain the same
        %                    as a CDF global attribute but
        %                    (1) it could be missing, or
        %                    (2) sometimes one may want to read an ROC-SGSE
        %                        dataset as if it was an RODP dataset or the other
        %                        way around.
        %
        function [OutputDatasetsMap] = produce_L2_LFR(InputDatasetsMap, rctDir, NsoTable, inputSciDsi, outputDsi, SETTINGS, L)
            
            InputHkCdf  = InputDatasetsMap('HK_cdf');
            InputCurCdf = InputDatasetsMap('CUR_cdf');
            InputSciCdf = InputDatasetsMap('SCI_cdf');

            %==============================
            % Configure bicas.calib object
            %==============================
            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inputSciDsi);
            useCtRcts = SETTINGS.get_fv('PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS')   && C.isL1r;
            useCti2   = SETTINGS.get_fv('PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2') && C.isL1r;
            
            if useCtRcts
                RctDataMap = bicas.calib.find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE(...
                    rctDir, 'LFR', ...
                    InputSciCdf.Ga.CALIBRATION_TABLE, ...
                    InputSciCdf.Zv.CALIBRATION_TABLE_INDEX, ...
                    InputSciCdf.Zv.BW, ...
                    L);
            else
                RctDataMap = bicas.calib.find_read_non_BIAS_RCTs_by_regexp(...
                    rctDir, SETTINGS, L);
            end
            Cal = bicas.calib(RctDataMap, rctDir, useCtRcts, useCti2, SETTINGS, L);
            
            %==============
            % Process data
            %==============
            HkSciTimePd           = bicas.proc_sub.process_HK_CDF_to_HK_on_SCI_TIME(InputSciCdf, InputHkCdf,  SETTINGS, L);
            InputSciCdf           = bicas.proc_sub.process_LFR_CDF_normalize(       InputSciCdf, inputSciDsi, SETTINGS, L);
            SciPreDc              = bicas.proc_sub.process_LFR_CDF_to_PreDC(        InputSciCdf, inputSciDsi, HkSciTimePd, SETTINGS, L);
            SciPostDc             = bicas.proc_sub.process_calibrate_demux(         SciPreDc, InputCurCdf, Cal,    SETTINGS, L);
            [SciPreDc, SciPostDc] = bicas.proc_sub.process_quality_filter(          SciPreDc, SciPostDc, NsoTable, SETTINGS, L);
            OutputSciCdf          = bicas.proc_sub.process_PostDC_to_LFR_CDF(       SciPreDc, SciPostDc, outputDsi, L);
            
            
            
            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('SCI_cdf') = OutputSciCdf;
        end



        % ARGUMENTS
        % =========
        % InputDatasetsMap : containers.Map: key=<argument key> --> value=PDV for input CDF
        %
        function [OutputDatasetsMap] = produce_L2_TDS(InputDatasetsMap, rctDir, NsoTable, inputSciDsi, outputDsi, SETTINGS, L)
            
            InputHkCdf  = InputDatasetsMap('HK_cdf');
            InputCurCdf = InputDatasetsMap('CUR_cdf');
            InputSciCdf = InputDatasetsMap('SCI_cdf');
            
            %==============================
            % Configure bicas.calib object
            %==============================
            % NOTE: TDS L1R never uses CALIBRATION_TABLE_INDEX2
            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inputSciDsi);
            if C.isTdsCwf
                settingUseCt = 'PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS';
                rctTypeId    = 'TDS-CWF';
            else
                settingUseCt = 'PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS';
                rctTypeId    = 'TDS-RSWF';
            end
            useCtRcts = SETTINGS.get_fv(settingUseCt) && C.isL1r;
            useCti2   = false;    % Always false for TDS.
            
            if useCtRcts
                RctDataMap = bicas.calib.find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE(...
                    rctDir, rctTypeId, ...
                    InputSciCdf.Ga.CALIBRATION_TABLE, ...
                    InputSciCdf.Zv.CALIBRATION_TABLE_INDEX, ...
                    [], ...
                    L);
            else
                RctDataMap = bicas.calib.find_read_non_BIAS_RCTs_by_regexp(...
                    rctDir, SETTINGS, L);
            end
            Cal = bicas.calib(RctDataMap, rctDir, useCtRcts, useCti2, SETTINGS, L);
            
            %==============
            % Process data
            %==============
            HkSciTimePd           = bicas.proc_sub.process_HK_CDF_to_HK_on_SCI_TIME(InputSciCdf, InputHkCdf,  SETTINGS, L);
            InputSciCdf           = bicas.proc_sub.process_TDS_CDF_normalize(       InputSciCdf, inputSciDsi, SETTINGS, L);
            SciPreDc              = bicas.proc_sub.process_TDS_CDF_to_PreDC(        InputSciCdf, inputSciDsi, HkSciTimePd, SETTINGS, L);
            SciPostDc             = bicas.proc_sub.process_calibrate_demux(         SciPreDc, InputCurCdf, Cal, SETTINGS, L);
            [SciPreDc, SciPostDc] = bicas.proc_sub.process_quality_filter(          SciPreDc, SciPostDc, NsoTable, SETTINGS, L);
            OutputSciCdf          = bicas.proc_sub.process_PostDC_to_TDS_CDF(       SciPreDc, SciPostDc, outputDsi, L);

            
            
            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('SCI_cdf') = OutputSciCdf;
        end
        
        
        
        function [OutputDatasetsMap] = produce_L3(InputDatasetsMap, NsoTable, SETTINGS, L)
            % Always the same DATASET_ID.
            INPUT_DATASET_ID = 'SOLO_L2_RPW-LFR-SURV-CWF-E';

            InputLfrCwfCdf = InputDatasetsMap('LFR-SURV-CWF-E_cdf');

            [InputLfrCwfCdf.Zv, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(InputLfrCwfCdf.Zv, ...
                {{{'VDC', 'V'}, 'VDC'}}, 'Assert one matching candidate');
            
            bicas.proc_sub.handle_zv_name_change(...
                fnChangeList, INPUT_DATASET_ID, SETTINGS, L, 'VDC', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')
            
            %===================================================================
            % Calculate
            % (1) E-field, and
            % (2) s/c potentials
            % via BICAS-external code
            % -----------------------
            % NOTE: Needs to be careful with the units, and incompatible updates
            % to solo.vdccal without the knowledge of the BICAS author.
            % Therefore extra assertions to detect such changes.
            %===================================================================
            TsVdc = TSeries(...
                EpochTT(InputLfrCwfCdf.Zv.Epoch), InputLfrCwfCdf.Zv.VDC, ...
                'TensorOrder', 1, ...
                'repres', {'x', 'y', 'z'});
            [TsEdc, TsPsp, TsScpot] = solo.vdccal(TsVdc);
            EJ_library.assert.sizes(...
                InputLfrCwfCdf.Zv.Epoch, [-1, 1], ...
                TsEdc.data,   [-1, 3], ...
                TsPsp.data,   [-1, 1], ...
                TsScpot.data, [-1, 3])
            assert(strcmp(TsEdc.units,   'mV/m'))
            assert(strcmp(TsPsp.units,   'V'))
            assert(strcmp(TsScpot.units, 'V'))
            
            %=================
            % Convert E-field
            %=================
            zvEdcMvpm = TsEdc.data;    % MVPM = mV/m
            % Set E_x = NaN, but only if assertion deems that the corresponding
            % information is missing
            % -----------------------------------------------------------------
            % IMPLEMENTATION NOTE: solo.vdccal set antenna 1 to be zero, if the
            % source data is non-fill value/NaN, but NaN if fill value. Must
            % therefore check for both zero and NaN.
            % Ex: 2020-08-01
            % IMPLEMENTATION NOTE: ismember does not work for NaN.
            assert(all(zvEdcMvpm(:, 1) == 0 | isnan(zvEdcMvpm(:, 1))), ...
                ['EDC for antenna 1 returned from BICAS_external code', ...
                ' solo.vdccal() is not zero or NaN and can therefore not be', ...
                ' assumed to unknown anymore. BICAS needs to be updated to reflect this.'])
            zvEdcMvpm(:, 1) = NaN;
            clear TsEdc
            
            
            
            EfieldCdf = struct();
            EfieldCdf.Epoch            = InputLfrCwfCdf.Zv.Epoch;
            EfieldCdf.QUALITY_BITMASK  = InputLfrCwfCdf.Zv.QUALITY_BITMASK;
            EfieldCdf.QUALITY_FLAG     = min(...
                InputLfrCwfCdf.Zv.QUALITY_FLAG, ...
                SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX'), ...
                'includeNaN');
            EfieldCdf.DELTA_PLUS_MINUS = InputLfrCwfCdf.Zv.DELTA_PLUS_MINUS;
            EfieldCdf.EDC_SFR          = zvEdcMvpm;
            
            
            
            ScpotCdf = struct();
            ScpotCdf.Epoch             = InputLfrCwfCdf.Zv.Epoch;
            ScpotCdf.QUALITY_BITMASK   = InputLfrCwfCdf.Zv.QUALITY_BITMASK;
            ScpotCdf.QUALITY_FLAG      = min(...
                InputLfrCwfCdf.Zv.QUALITY_FLAG, ...
                SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX'), ...
                'includeNaN');
            ScpotCdf.DELTA_PLUS_MINUS  = InputLfrCwfCdf.Zv.DELTA_PLUS_MINUS;
            ScpotCdf.SCPOT             = TsScpot.data;
            ScpotCdf.PSP               = TsPsp.data;

            
            
            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('EFIELD_cdf') = EfieldCdf;
            OutputDatasetsMap('SCPOT_cdf')  = ScpotCdf;
        end
        
        
        
    end    % methods(Static, Access=public)
    
    
    
end    % classdef
