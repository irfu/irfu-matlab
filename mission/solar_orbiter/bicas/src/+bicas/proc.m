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
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
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
    % 
    % TODO-DECISION: Use PDID system?
    %   NOTE: data_manager_old's PDID used skeleton versions.
    %   NOTE: According to RCS ICD 00037 iss1/rev2, draft 2019-07-11, the s/w descriptor interface no longer specifies
    %         the version of the input datasets. One can still specify modes (and CLI parameters) that require specific
    %         input skeleton versions though.
    %   NOTE: Processing functions still need to know the input skeleton version.
    %   NOTE: One can choose PDIDs such that the incorporate the version or not. If they do not, then the corresponding
    %         PDVs must themselves contain the same information.
    % --
    % PROPOSAL: Production functions should not assume/specify any particular input dataset version, but read it out
    %           from global attributes (part of the PDV).
    % PROPOSAL: Somehow associate metadata with each function.
    %   PRO: For a given OUTPUT dataset/PDID, one can get the list of possible input datasets/PDIDs
    %   PROPOSAL: A given production function with a particular output dataset/PDID represents a s/w mode and will (for
    %       the right arguments, e.g. empty) return information on the corresponding inputs (incl. for s/w descriptor).
    % 
    %   NOTE: Metadata associated with each s/w mode
    %
    % --
    % TODO-DECISION: How handle differences between pipelines?
    %
    % PROPOSAL: (Outside) Derive PDID from input dataset, then use it when sending input dataset argument to production function.
    %   CON: Possible to derive non-existant (non-supported) PDIDs. Not the best way to test input datasets.

    methods(Static, Access=public)
        
        % ARGUMENTS
        % =========
        % InputDatasetsMap : containers.Map: key=<argument key> --> value=PDV for input CDF
        % inputSciDsi      : The science input dataset will be interpreted as having this DATASET_ID.
        %                    RATIONALE: InputDatasetsMap should contain the same as a CDF global attribute but
        %                    (1) it could be missing, or
        %                    (2) sometimes one may want to read an ROC-SGSE dataset as if it was an RODP dataset or the other way around.
        %
        function [OutputDatasetsMap] = produce_L2_LFR(InputDatasetsMap, Cal, inputSciDsi, outputDsi, outputVersion, SETTINGS, L)
            
            InputHkPd  = InputDatasetsMap('HK_cdf');
            InputCurPd = InputDatasetsMap('CUR_cdf');
            InputSciPd = InputDatasetsMap('SCI_cdf');

            %==============================
            % Configure calibration object
            %==============================
            C = EJ_library.so.adm.classify_DATASET_ID(inputSciDsi);
            useCt   = SETTINGS.get_fv('PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS')   && C.isL1R;
            useCti2 = SETTINGS.get_fv('PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2') && C.isL1R;
            if useCt
                Cal.read_non_BIAS_RCT_by_CALIBRATION_TABLE('LFR', ...
                    InputSciPd.Ga.CALIBRATION_TABLE, ...
                    InputSciPd.Zv.CALIBRATION_TABLE_INDEX, ...
                    useCti2);
            else
                Cal.read_non_BIAS_RCTs_by_regexp(useCti2);
            end
            
            HkSciTimePd  = bicas.proc_sub.process_HK_to_HK_on_SCI_TIME(  InputSciPd, InputHkPd,  SETTINGS, L);
            %CurSciTimePd = bicas.proc_sub.process_CUR_to_CUR_on_SCI_TIME(InputSciPd, InputCurPd, SETTINGS, L);
            SciPreDcPd   = bicas.proc_sub.process_LFR_to_PreDC(          InputSciPd, inputSciDsi, HkSciTimePd, SETTINGS, L);
            SciPostDcPd  = bicas.proc_sub.process_demuxing_calibration(  SciPreDcPd, InputCurPd, Cal, SETTINGS, L);
            OutputSciPd  = bicas.proc_sub.process_PostDC_to_LFR(         SciPostDcPd, outputDsi);
            
            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('SCI_cdf') = OutputSciPd;
        end



        % ARGUMENTS
        % =========
        % InputDatasetsMap : containers.Map: key=<argument key> --> value=PDV for input CDF
        %
        function [OutputDatasetsMap] = produce_L2_TDS(InputDatasetsMap, Cal, inputSciDsi, outputDsi, outputVersion, SETTINGS, L)
            
            InputHkPd  = InputDatasetsMap('HK_cdf');
            InputCurPd = InputDatasetsMap('CUR_cdf');
            InputSciPd = InputDatasetsMap('SCI_cdf');
            
            %==============================
            % Configure calibration object
            %==============================
            % NOTE: TDS L1R never uses CALIBRATION_TABLE_INDEX2
            C = EJ_library.so.adm.classify_DATASET_ID(inputSciDsi);
            if C.isTdsCwf
                settingUseCt   = 'PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS';
                %settingUseCti2 = 'PROCESSING.L1R.TDS.CWF.USE_ZV_CALIBRATION_TABLE_INDEX2';
                rctId          = 'TDS-CWF';
            else
                settingUseCt   = 'PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS';
                %settingUseCti2 = 'PROCESSING.L1R.TDS.RSWF.USE_ZV_CALIBRATION_TABLE_INDEX2';
                rctId          = 'TDS-RSWF';
            end
            useCt   = SETTINGS.get_fv(settingUseCt)   && C.isL1R;
            %useCti2 = SETTINGS.get_fv(settingUseCti2) && C.isL1R;
            useCti2 = false;
            if useCt
                Cal.read_non_BIAS_RCT_by_CALIBRATION_TABLE(rctId, ...
                    InputSciPd.Ga.CALIBRATION_TABLE, ...
                    InputSciPd.Zv.CALIBRATION_TABLE_INDEX, ...
                    useCti2);
            else
                Cal.read_non_BIAS_RCTs_by_regexp(useCti2);
            end
            
            
            
            HkSciTimePd = bicas.proc_sub.process_HK_to_HK_on_SCI_TIME(InputSciPd, InputHkPd, SETTINGS, L);
            SciPreDcPd  = bicas.proc_sub.process_TDS_to_PreDC(        InputSciPd, inputSciDsi, HkSciTimePd, SETTINGS, L);
            SciPostDcPd = bicas.proc_sub.process_demuxing_calibration(SciPreDcPd, InputCurPd, Cal, SETTINGS, L);
            OutputSciPd = bicas.proc_sub.process_PostDC_to_TDS(       SciPostDcPd, outputDsi);

            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('SCI_cdf') = OutputSciPd;

        end
        
    end    % methods(Static, Access=public)
    
end    % classdef
