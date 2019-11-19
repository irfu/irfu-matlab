%
% Singleton class that stores (after having "built" tit) an unmodifiable data structure that represents which and how
% s/w modes are CURRENTLY VISIBLE to the user. What that data structure contains thus depends on
% -- current pipeline: RODP, ROC-SGSE
% -- whether support for old L2R input datasets is enabled or not.
%
% Data here
% -- is intended to (1) add metadata to the production functions for
%       (a) the user interface, and
%       (b) the s/w descriptor
% -- contains variables to make it possible to match information for input and output datasets here, with that of
%    bicas.proc' production functions.
%
%
%
% IMPLEMENTATION NOTES
% ====================
% The class essentially consists of one large struct, and a constructor that builds it. The large data struct contains
% many parts which are similar but not the same. To do this, much of the data is "generated" with hardcoded strings
% (mostly the same in every iteration), in which specific codes/substrings are substituted algorithmically (different in
% different iterations). To avoid mistakes, the code uses a lot of assertions to protect against mistakes, e.g.
% -- algorithmic bugs
% -- mistyped hardcoded info
% -- mistakenly confused arguments with each other.
% Assertions are located at the place where "values are placed in their final location".
% 
% NOTE: To implement backward compatibility with L2R input datasets, the code must be able to handle
% -- changing input dataset levels: L1R (new), L2R (old).
% -- DATASET_ID with (new) and without (old) a trailing "-E".
% It implements L2R input datasets via separate S/W modes.
% 
% RATIONALE:
% -- Should decrease the amount of overlapping hardcoded information to e.g. reduce risk of mistakes, reduce manual
%    work when verifying updates.
% -- Having one big, somewhat redundant data structure should make the interface to the rest of BICAS relatively
%    future-proof, in the face of future updates
% -- Useful for expected future bias current datasets
% -- Possible need for backward compatibility
%
%
% DEFINITIONS
% ===========
% SIP = "Specific Input Parameters" (RCS ICD).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-31
%
classdef swmode_defs
    % PROPOSAL: New class name implying that it only contains S/W modes VISIBLE to the user, that it DEFINES what is
    % visible.
    %
    % PROPOSAL: Pick SWD name/descriptions from master CDFs.
    % PROPOSAL: Obtain output dataset level from production function metadata?!!
    % PROPOSAL: Include output dataset version.
    %   PRO: Can deduce name of master CDF.
    %       PRO: Needed for SWD.
    %       PRO: Needed for deriving master-CDF filename.
    %   PRO: Needed for verifying conformance with production function.
    %
    % PROPOSAL: Always produce all possible s/w modes (both pipelines, incl. L2R), then filter out the undesired ones
    % using internal metadata for every S/W mode.
    %
    % PROPOSAL: Use PF = prodFunc, production function
    % PROPOSAL: Same input CDF can have multiple DATASET_IDs, but only one is shown in the s/w descriptor.
    %   PRO: Can handle old datasets with ROG-SGSE DATASET_IDs, and otherwise only use RODP DATASET_IDs.



    % PRIVATE, STATIC, CONSTANTS
    properties(Constant, GetAccess=private)
        % The RCS ICD 00037, iss1rev2, draft 2019-07-11, section 5.3 seems (ambiguous) to imply this regex for S/W mode CLI parameters.
        % regexp: "\w    A word character [a-z_A-Z0-9]"
        SW_MODE_CLI_PARAMETER_REGEX = '^[A-Za-z][\w-]+$';   % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.
        
        % The RCS ICD 00037 iss1rev2 draft 2019-07-11, section 3.1.2.3 only permits these characters (and only lowercase!).
        % This regexp only describes the "option body", i.e. not the preceding "--".
        SIP_CLI_OPTION_BODY_REGEX = '[a-z0-9_]+';
    end

    
    
    % PUBLIC, IMMUTABLE
    properties(SetAccess=immutable)
        List
    end
    
    
    
    % PRIVATE, INSTANCE, IMMUTABLE
    properties(GetAccess=private, SetAccess=immutable)
        dsiPipelinePrefix    % Prefix in DATASET_ID (DSI).
        outputDatasetLevel
    end



    methods(Access=public)
        
        % Constructor
        %
        % ARGUMENTS
        % =========
        % enableRocsgseL2rInput : true/false, 1/0. Whether to enable (make visible) support for ROC-SGSE.
        %
        % IMPLEMENTATION NOTE: The constructor used to be written so that it was easy to disable S/W modes with L2R
        % input datasets (for backward compatibility). That functionality has now been now removed, although the
        % implementation has not been entirely updated to take advantage of this (not simplified of this).
        % 
%        function obj = swmode_defs(pipelineId, enableRocsgseL2rInput, enableTds)
        function obj = swmode_defs(pipelineId, enableTds)
            % PROPOSAL: Re-implement (top-level) hard-coded constants by setting multiple redundant 1D(?) vectors that covers every case.
            %   Then set various cases by assigning constants to many elements using MATLAB syntax.
            %   One index representing: Combination of DATASET_ID+Skeleton_Version (both pipelines, LFR+TDS, HK+SCI), every element contains data for
            %   that dataset. Must use combination DATASET_ID+Skeleton_Version to potentially cover old versions.            
            %   Manipulate and set multiple elements smoothly by using vectors for indices.
            %   Ex: Vectors to set: pipelineIdVector, skeletonVersionVector, SBMx_SURV_vector, CWF_SWF_vector, output
            %       dataset level, vectors for human-readable description string(s) (e.g. modeStr)
            %   Ex: Vectors with indices for dataset in/for: either pipeline, science or HK, LFR or TDS, latest versions
            %       or backward-compatibility versions.
            %   --
            %   TODO-DECISION: Above describes input & output (?) data sets. How relates to s/w modes?
            %   
            
            %========================================
            % Select constants depending on pipeline
            %========================================
            switch(pipelineId)
                case {'ROC-SGSE', 'RGTS'}
                    %obj.dsiPipelinePrefix     = 'ROC-SGSE';      % Prefix in DATASET_ID (DSI).
                    obj.dsiPipelinePrefix     = 'SOLO';       % Prefix in DATASET_ID (DSI). TEST: On request from ROC.
                    inputDatasetLevelList     = {'L1R'};
%                    obj.outputDatasetLevel    =  'L2S';
                    obj.outputDatasetLevel    =  'L2';    % On request from ROC.
                    inputDashEList            = {'-E'};
                    swModeCliOptionAmendmList = {''};
                    lfrOutputSkeletonVersion  = {'05', '05', '05', '05'};
                    tdsOutputSkeletonVersion  = {'05', '05'};
                    
%                     if enableRocsgseL2rInput
%                         inputDatasetLevelList{end+1}     = 'L2R';     % NOTE: L2R etc only kept for backward-compatibility.
%                         inputDashEList{end+1}            = '';
%                         swModeCliOptionAmendmList{end+1} = '_L2R';
%                     end
                    
                case 'RODP'
                    obj.dsiPipelinePrefix     = 'SOLO';    % NOTE: SOLO, not RODP.
                    inputDatasetLevelList     = {'L1R'};
                    inputDashEList            = {'-E'};
                    swModeCliOptionAmendmList = {''};
                    obj.outputDatasetLevel    = 'L2';
                    lfrOutputSkeletonVersion  = {'05', '05', '05', '05'};
                    tdsOutputSkeletonVersion  = {'05', '05'};

                otherwise
                    error('BICAS:swmode_defs:Assertion:IllegalArgument', 'Can not interpret "pipelineId=%s', pipelineId)
            end
            clear pipelineId
            

            
            % Define function which interprets (replaces) specific substrings.            
            % "strmod" = string modify, "g"=global
            strmodg = @(s, iInputLevel) strrep(strrep(strrep(strrep(strrep(s, ...
                '<PLP>', obj.dsiPipelinePrefix), ...
                '<LI>',  inputDatasetLevelList{iInputLevel}), ...
                '<LO>',  obj.outputDatasetLevel), ...
                '<I-E>', inputDashEList{iInputLevel}), ...
                '<L2R amendm>', swModeCliOptionAmendmList{iInputLevel});    % I-E: "I"=Input, "-E"=optional "-E"
            
            % Input def that is reused multiple times.
            HK_INPUT_DEF = obj.def_input_dataset(...
                'in_hk', strmodg('<PLP>_HK_RPW-BIA', 1), 'HK_cdf');



            LFR_SW_MODE_DATA = struct(...
                'SBMx_SURV',       {'SBM1', 'SBM2', 'SURV', 'SURV'}, ...
                'CWF_SWF',         {'CWF',  'CWF',  'CWF',  'SWF'}, ...
                'modeStr',         {'selective burst mode 1', 'selective burst mode 2', 'survey mode', 'survey mode'}, ...
                'outputSkeletonVersion', lfrOutputSkeletonVersion);
            TDS_SW_MODE_DATA = struct(...
                'CWF_RSWF',        {'CWF', 'RSWF'}, ...
                'outputSkeletonVersion', tdsOutputSkeletonVersion);
            
            
            
            List = struct('prodFunc', {}, 'cliOption', {}, 'swdPurpose', {}, 'inputsList', {}, 'outputsList', {});
            for iInputLevel = 1:numel(inputDatasetLevelList)
                
                %==============================================
                % Iterate over the "fundamental" LFR S/W modes
                %==============================================
                for iSwm = 1:length(LFR_SW_MODE_DATA)
                    strmod = @(s) strrep(strrep(strrep(strmodg(s, iInputLevel), ...
                        '<SBMx/SURV>',  LFR_SW_MODE_DATA(iSwm).SBMx_SURV), ...
                        '<C/SWF>',      LFR_SW_MODE_DATA(iSwm).CWF_SWF), ...
                        '<mode str>',   LFR_SW_MODE_DATA(iSwm).modeStr);
                    
                    SCI_INPUT_DEF = obj.def_input_dataset(...
                        'in_sci', ...
                        strmod('<PLP>_<LI>_RPW-LFR-<SBMx/SURV>-<C/SWF><I-E>'), ...
                        'SCI_cdf');
                    
                    SCI_OUTPUT_DEF = obj.def_output_dataset(...
                        strmod('<PLP>_<LO>_RPW-LFR-<SBMx/SURV>-<C/SWF>-E'), ...
                        strmod('LFR <LO> <C/SWF> science electric <mode str> data'), ...
                        strmod('RPW LFR <LO> <C/SWF> science electric (potential difference) data in <mode str>, time-tagged'), ...
                        LFR_SW_MODE_DATA(iSwm).outputSkeletonVersion);

                    List(end+1) = obj.def_swmode(...
                        @(InputDatasetsMap, Cal) bicas.proc.produce_L2_LFR(...
                            InputDatasetsMap, ...
                            Cal, ...
                            SCI_INPUT_DEF.datasetId, ...
                            SCI_OUTPUT_DEF.datasetId, ...
                            SCI_OUTPUT_DEF.skeletonVersion), ...
                        strmod('LFR-<SBMx/SURV>-<C/SWF>-E<L2R amendm>'), ...
                        strmod('Generate <SBMx/SURV> <C/SWF> electric field <LO> data (potential difference) from LFR <LI> data'), ...
                        [SCI_INPUT_DEF, HK_INPUT_DEF], [SCI_OUTPUT_DEF]);
                end
                
                if enableTds
                    %==============================================
                    % Iterate over the "fundamental" TDS S/W modes
                    %==============================================
                    for iSwm = 1:numel(TDS_SW_MODE_DATA)
                        strmod = @(s) strrep(strmodg(s, iInputLevel), ...
                            '<C/RSWF>', TDS_SW_MODE_DATA(iSwm).CWF_RSWF);
                        
                        SCI_INPUT_DEF = obj.def_input_dataset(...
                            'in_sci', ...
                            strmod('<PLP>_<LI>_RPW-TDS-LFM-<C/RSWF><I-E>'), ...
                            'SCI_cdf');
                        
                        SCI_OUTPUT_DEF = obj.def_output_dataset(...
                            strmod('<PLP>_<LO>_RPW-TDS-LFM-<C/RSWF>-E'), ...
                            strmod('LFR <LO> <C/RSWF> science electric LF mode data'), ...
                            strmod('RPW LFR <LO> <C/RSWF> science electric (potential difference) data in LF mode, time-tagged'), ...
                            TDS_SW_MODE_DATA(iSwm).outputSkeletonVersion);
                        
                        List(end+1) = obj.def_swmode(...
                            @(InputDatasetsMap, Cal) bicas.proc.produce_L2_TDS(...
                                InputDatasetsMap, ...
                                Cal, ...
                                SCI_INPUT_DEF.datasetId, ...
                                SCI_OUTPUT_DEF.datasetId, ...
                                SCI_OUTPUT_DEF.skeletonVersion), ...
                            strmod('TDS-LFM-<C/RSWF>-E<L2R amendm>'), ...
                            strmod('Generate <C/RSWF> electric field <LO> data (potential difference) from TDS LF mode <LI> data'), ...
                            [SCI_INPUT_DEF, HK_INPUT_DEF], [SCI_OUTPUT_DEF]);
                    end
                end
            end    % for iInputLevel = 1:numel(inputDatasetLevelList)
            
            
            
            obj.List = List;
            clear List
            
            EJ_library.utils.assert.castring_set({obj.List(:).cliOption})
        end    % Constructor



        function swModeInfo = get_sw_mode_info(obj, swModeCliOption)
            i = find(strcmp(swModeCliOption, {obj.List(:).cliOption}));
            EJ_library.utils.assert.scalar(i)
            swModeInfo = obj.List(i);
        end
        
    end    % methods(Access=public)
    
    

    methods(Access=private)
        
        % NOTE: Could technically be a static method. Only instance method for grouping it with other analogous methods.
        function Def = def_swmode(~, prodFunc, cliOption, swdPurpose, inputsList, outputsList)
            Def.prodFunc    = prodFunc;
            Def.cliOption   = cliOption;   % NOTE: s/w mode CLI _ARGUMENT_ is not intended to be prefixed by e.g. "--". Variable therefore NOT named *Body.
            Def.swdPurpose  = swdPurpose;
            Def.inputsList  = inputsList;
            Def.outputsList = outputsList;
            
            
            
            % ASSERTIONS
            EJ_library.utils.assert.castring_set( {...
                Def.inputsList(:).cliOptionHeaderBody, ...
                Def.outputsList(:).cliOptionHeaderBody })   % Important. Check uniqueness of SIP options.
            EJ_library.utils.assert.castring_set( {...
                Def.inputsList(:).prodFuncInputKey })   % Maybe not really necessary.
            EJ_library.utils.assert.castring_set( {...
                Def.outputsList(:).prodFuncOutputKey })   % Maybe not really necessary.
            
            bicas.swmode_defs.assert_SW_mode_CLI_option(Def.cliOption)
            bicas.swmode_defs.assert_text(              Def.swdPurpose)
            
            EJ_library.utils.assert.castring_set({...
                Def.inputsList(:).cliOptionHeaderBody, ...
                Def.outputsList(:).cliOptionHeaderBody})
        end

        
        
        function Def = def_input_dataset(obj, cliOptionHeaderBody, datasetId, prodFuncInputKey)
            % NOTE: No dataset/skeleton version.
            Def.cliOptionHeaderBody = cliOptionHeaderBody;
            Def.prodFuncInputKey    = prodFuncInputKey;
            Def.datasetId           = datasetId;
            
            bicas.swmode_defs.assert_SIP_CLI_option(Def.cliOptionHeaderBody)
            obj.assert_DATASET_ID(                  Def.datasetId)    % NOTE: Using the internal assertion function, not the global one.
        end

        
        
        function Def = def_output_dataset(obj, datasetId, swdName, swdDescription, skeletonVersion)
            Def.cliOptionHeaderBody = 'out_sci';
            Def.prodFuncOutputKey   = 'SCI_cdf';
            Def.swdName             = swdName;
            Def.swdDescription      = swdDescription;
            Def.datasetId           = datasetId;
            Def.datasetLevel        = obj.outputDatasetLevel;     % NOTE: Automatically set.
            Def.skeletonVersion     = skeletonVersion;
            
            bicas.swmode_defs.assert_SW_mode_CLI_option(Def.cliOptionHeaderBody)
            bicas.swmode_defs.assert_text(              Def.swdName)
            bicas.swmode_defs.assert_text(              Def.swdDescription)
            obj.assert_DATASET_ID(                      Def.datasetId)
            bicas.assert_dataset_level(                 Def.datasetLevel)
            bicas.assert_skeleton_version(              Def.skeletonVersion)
        end



        % NOTE: Wrapper around global counterpart.
        function assert_DATASET_ID(obj, datasetId)
            bicas.assert_DATASET_ID(datasetId)
            
            % ASSERTION: Pipeline
            assert(strcmp(...
                obj.dsiPipelinePrefix, ...
                datasetId(1:numel(obj.dsiPipelinePrefix))...
                ))
        end

    end    % methods(Access=private)    



    methods(Static, Access=private)

        % Assert that string contains human-readable text.
        function assert_text(str)
            EJ_library.utils.assert.castring_regexp(str, '.* .*')
            EJ_library.utils.assert.castring_regexp(str, '[^<>]*')
        end
        
        function assert_SW_mode_CLI_option(swModeCliOption)
            EJ_library.utils.assert.castring_regexp(swModeCliOption, bicas.swmode_defs.SW_MODE_CLI_PARAMETER_REGEX)
        end

        % NOTE: Really refers to "option body".
        function assert_SIP_CLI_option(sipCliOptionBody)
            EJ_library.utils.assert.castring_regexp(sipCliOptionBody, bicas.swmode_defs.SIP_CLI_OPTION_BODY_REGEX)
        end
        
    end    % methods(Static, Access=private)

end
