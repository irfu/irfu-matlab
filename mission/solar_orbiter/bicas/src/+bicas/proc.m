%
% Set of functions for producing one specific output dataset PDV from the necessary input dataset PDVs.
%
%
% PRODUCTION FUNCTIONS
% ====================
% A function with interface
%   OutputsMap = produce_*(InputsMap)
% where
% InputsMap  : containers.Map with
%                <keys>       : String defining a name of an input ("prodFuncInputKey" in swmode_defs).
%                <values>     : A struct with data corresponding to a CDF file (zVariables+global attributes).
% OutputsMap : containers.Map with
%                <keys>       : String defining a name of an output ("prodFuncOutputKey" in swmode_defs).
%                <values>     : A struct with data corresponding to a CDF file (zVariables).
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
    %
    % TODO-NEED-INFO: Do the production functions really need to be aware of EIn & EOut PDIDs?!
    %       They need to know which input is which though: HK vs SCI, different types of datasets.
    % --
    % PROPOSAL: Submit CDF global attributes to processing functions.
    %     PRO: Can use the ~Skeleton_version GA for assertion & interpreting data instead of PDID.
    %         CON: ~Skeleton_version GA can be wrong.
    %             Ex: Global attribute Skeleton_version
    %             PROPOSAL: Setting for overriding global attribute dataset version.
    %                 NOTE: Such setting needs two variables in principle:
    %                     --DATASET_ID for which Dataset_version will be overwritten,
    %                     --Dataset_version itself
    %                 PROPOSAL: Incorporate DATASET_ID in settings key.
    %                     Ex: ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_Dataset_version_OVERRIDE
    %                     CON: Presupposes that such settings pre-exist for every DATASET_ID.
    %                 PROPOSAL: Use settings value arrays somehow. May need to implement.
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
    % PROPOSAL: Master CDF as input.
    %   PRO: Needed for merging global attributes?
    %
    % PROPOSAL: Instantiate class, use instance methods instead of static.
    %   PRO: Can have SETTINGS & constants as instance variable instead of calling global variables.
    %
    % TODO-DECISION: How handle differences between pipelines?
    %
    % PROPOSAL: (Outside) Derive PDID from input dataset, then use it when sending input dataset argument to production function.
    %   CON: Possible to derive non-exististant (non-supported) PDIDs. Not the best way to test input datasets.
    
    methods(Static, Access=public)
        
        % Produce PD for PDID=V03_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E.
        %
        % ARGUMENTS
        % =========
        % InputsMap : containers.Map: key=<argument key> --> value=PDV for input CDF
        %
        function [OutputsMap] = produce_L2S_L2_LFR(InputsMap, outputDsi, outputVersion)
               
                HkPd  = InputsMap('HK_cdf');
                SciPd = InputsMap('SCI_cdf');
                
                HkSciTimePd = bicas.proc_sub.process_HK_to_HK_on_SCI_TIME(SciPd, HkPd);
                
                SciPreDcPd  = bicas.proc_sub.process_LFR_to_PreDC(        SciPd, HkSciTimePd);
                SciPostDcPd = bicas.proc_sub.process_demuxing_calibration(SciPreDcPd);
                OutputSciPd = bicas.proc_sub.process_PostDC_to_LFR(       SciPostDcPd, outputDsi, outputVersion);
                
                OutputsMap = containers.Map();
                OutputsMap('SCI_cdf') = OutputSciPd;
%             end
        end
        
    end    % methods(Static, Access=public)
    
end    % classdef
