%
% Set of functions for producing one specific output dataset PDV from the necessary input dataset PDVs.
%
%
% PRODUCTION FUNCTIONS
% ====================
% A function with interface
%   ProcessData = produce_*(InputsMap)
% where
% InputsMap : containers.Map with
%             InputsMap : A containers.Map
%                <keys>       : PDID.
%                <values>     : A struct with fields .pd (process data) and .pdid (PDID for .pd).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-30
%
classdef pipelines
    % PROPOSAL: Other name of class.
    %   PROPOSAL: production_functions
    %   PROPOSAL: production
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
    %       CON: Conflicts with RODP, ROC-SGSE pipelines.
    %   PROPOSAL: swmode_pipelines
    %   TODO-DECISION: Include skeleton version in function names?
    %   
    % TODO-DECISION: Use PDID system?
    %   NOTE: data_manager_old's PDID uses skeleton versions.
    %   NOTE: According to RCS ICD 00037 iss1/rev2, draft 2019-07-11, the s/w descriptor interface no specifies
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
    %             Ex: Setting INPUT_CDF_ASSERTIONS.STRICT_SKELETON_VERSION.
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
    %   PRO: Can have SETTINGS, CONSTANTS as instance variable instead of calling global variables.
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
        %function [OutputsMap, InputPdidsList, OutputPdidList] = produce_L2S_L2_LFR(InputsMap, outputPdid, inputLevel, pipelineId)
        function [OutputsMap] = produce_L2S_L2_LFR(InputsMap, outputDsi, outputVersion)
            
%             OutputPdidList = {outputPdid};
%             
%             switch(pipelineId)
%                 case {'ROC-SGSE', 'RGTS'}
%                     dsiPrefix = 'ROC-SGSE';
%                     switch(inputLevel)
%                         case 'L2R'
%                             dashESuffix = '';
%                         case 'L1R'
%                             dashESuffix = '-E';
%                         otherwise
%                             error('BICAS:pipelines:produce_L2S_L2_LFR:Assertion:IllegalArgument', 'Illegal inputLevel="%s"', inputLevel)
%                     end
%                 case 'RODP'
%                     dsiPrefix   = 'SOLO';
%                     dashESuffix = '-E';
%                     assert(strcmp(inputLevel, 'L1R'))
%                 otherwise
%                     error('BICAS:pipelines:produce_L2S_L2_LFR:Assertion:IllegalArgument', 'Illegal pipelineId="%s"', pipelineId)
%             end
%             
%             % BUG: outputPdid can only be ROC-SGSE.
%             switch(outputPdid)
%                 case 'V03_ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E'
%                     InputSciPdid = '<PLP>_<LI>_RPW-LFR-SBM1-CWF<-E>';
%                 case 'V03_ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E'
%                     InputSciPdid = '<PLP>_<LI>_RPW-LFR-SBM2-CWF<-E>';
%                 case 'V03_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E'
%                     InputSciPdid = '<PLP>_<LI>_RPW-LFR-SURV-CWF<-E>';
%                 case 'V03_ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E'
%                     InputSciPdid = '<PLP>_<LI>_RPW-LFR-SURV-SWF<-E>';
%                 otherwise
%                     error('pipelines:produce_L2S_L2_LFR:Assertion:IllegalCodeConfiguration', ...
%                         'This function can not handle outputPdid="%s"', outputPdid)
%             end
%             InputSciPdid = strrep(InputSciPdid, '<PLP>', dsiPrefix);
%             InputSciPdid = strrep(InputSciPdid, '<LI>',  inputLevel);
%             InputSciPdid = strrep(InputSciPdid, '<-E>',  dashESuffix);            
%             
%             InputHkPdid = '<PLP>_HK_RPW-BIA';
%             InputHkPdid = strrep(InputHkPdid, '<PLP>', dsiPrefix);
%             
%             InputPdidsList = {InputSciPdid, InputHkPdid};

%             if isempty(InputsMap)
%                 %===============================
%                 % CASE: No processing requested
%                 %===============================
%                 OutputsMap = [];
%             else
                %============================
                % CASE: Processing requested
                %============================
                
                HkPd  = InputsMap('HK_cdf');
                SciPd = InputsMap('SCI_cdf');
                
                HkSciTimePd = bicas.dm_processing_functions.process_HK_to_HK_on_SCI_TIME(SciPd, HkPd);
                
                SciPreDcPd  = bicas.dm_processing_functions.process_LFR_to_PreDC(        SciPd, HkSciTimePd);
                SciPostDcPd = bicas.dm_processing_functions.process_demuxing_calibration(SciPreDcPd);
                OutputSciPd = bicas.dm_processing_functions.process_PostDC_to_LFR(       SciPostDcPd, outputDsi, outputVersion);
                
                OutputsMap = containers.Map();
                OutputsMap('SCI_cdf') = OutputSciPd;
%             end
        end
        
    end    % methods(Static, Access=public)
    
end    % classdef
