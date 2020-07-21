%
% Singleton class that stores (after having "built" it) an unmodifiable data structure that represents which and how
% s/w modes are CURRENTLY VISIBLE to the user. What that data structure contains thus depends on
% -- current pipeline: RODP, ROC-SGSE
% -- whether support for L1 input datasets is enabled or not.
%
% Data here
% -- is intended to (1) add metadata to the production functions for
%       (a) the caller interface
%       (b) potentially future help text
%       (c) the s/w descriptor
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
% NOTE: To implement compatibility with L1 input datasets, the code must be able to handle
% -- changing input dataset levels: L1 (inofficial support), L1R (official support).
% It implements support for L1 input datasets via separate S/W modes.
% 
%
% RATIONALE
% =========
% -- Should decrease the amount of overlapping hardcoded information to e.g. reduce risk of mistakes, reduce manual
%    work when verifying updates.
% -- Having one big, somewhat redundant data structure should make the interface to the rest of BICAS relatively
%    future-proof, in the face of future updates.
% -- Useful for expected future bias current datasets.
% -- Possible need for backward compatibility.
%
%
% DEFINITIONS
% ===========
% SIP = "Specific Input Parameters" (RCS ICD).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-07-31
%
classdef swmode_defs
    % PROPOSAL: New class name implying that it only contains S/W modes VISIBLE to the user, that it DEFINES what is
    % visible for a given BICAS run.
    %
    % PROPOSAL: Pick SWD name/descriptions from master CDFs.
    % PROPOSAL: Obtain output dataset level from production function metadata?!!
    % PROPOSAL: Include output dataset version.
    %   PRO: Can deduce name of master CDF.
    %       PRO: Needed for SWD.
    %       PRO: Needed for deriving master-CDF filename.
    %   PRO: Needed for verifying conformance with production function.
    %
    % PROPOSAL: Always produce all possible s/w modes (both pipelines, incl. L1), then filter out the undesired ones
    % using internal metadata for every S/W mode.
    %
    % PROPOSAL: Use PF = prodFunc, production function
    % PROPOSAL: Same input CDF can have multiple DATASET_IDs, but only one is shown in the s/w descriptor.
    %   PRO: Can handle old datasets with ROG-SGSE DATASET_IDs, and otherwise only use RODP DATASET_IDs.
    %
    % PROPOSAL: Use classes instead of structs.
    %   PROPOSAL: Use sub-package for classes.
    %   PROPOSAL: Class for s/w mode definition.
    %       PROPOSAL: Name swmode_def
    %           CON: Too similar to swmode_defs
    %   PROPOSAL: Class for input datasets.
    %   PROPOSAL: Class for output datasets.
    %



    % PRIVATE, STATIC, CONSTANTS
    properties(Constant, GetAccess=private)
        % The RCS ICD 00037, iss1rev2, draft 2019-07-11, section 5.3 seems (ambiguous) to imply this regex for S/W mode
        % CLI parameters.
        % regexp: "\w    A word character [a-z_A-Z0-9]"
        %
        % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.
        SW_MODE_CLI_PARAMETER_REGEX = '^[A-Za-z][\w-]+$';   
        
        % The RCS ICD 00037 iss1rev2 draft 2019-07-11, section 3.1.2.3 only permits these characters (and only
        % lowercase!).
        % This regexp only describes the "option body", i.e. not the preceding "--".
        SIP_CLI_OPTION_BODY_REGEX = '[a-z0-9_]+';
    end

    
    
    % PUBLIC, IMMUTABLE
    properties(SetAccess=immutable)
        
        % NOTE: Implicit that it is a list of s/w modes (not in name). Note that it is a public property.
        List
    end
    
    
    
    methods(Access=public)
        
        
        
        % Constructor
        %
        % ARGUMENTS
        % =========
        %
        % IMPLEMENTATION NOTE: The constructor used to be written so that it was easy to disable S/W modes with L2R
        % input datasets (for backward compatibility). That functionality has now been now removed, although the
        % implementation has not been entirely updated to take advantage of this (not simplified of this).
        % 
        function obj = swmode_defs(SETTINGS, L)
            % PROPOSAL: Re-implement (top-level) hard-coded constants by setting multiple redundant 1D(?) vectors that
            %   covers every case. Then set various cases by assigning constants to many elements using MATLAB syntax.
            %   One index representing: Combination of DATASET_ID+Skeleton_Version (both pipelines, LFR+TDS, HK+SCI),
            %   every element contains data for that dataset. Must use combination DATASET_ID+Skeleton_Version to
            %   potentially cover old versions.
            %   Manipulate and set multiple elements smoothly by using vectors for indices.
            %   Ex: Vectors to set: skeletonVersionVector, SBMx_SURV_vector, CWF_SWF_vector, output
            %       dataset level, vectors for human-readable description string(s) (e.g. modeStr)
            %   Ex: Vectors with indices for dataset in/for: either pipeline, science or HK, LFR or TDS, latest versions
            %       or backward-compatibility versions.
            %   --
            %   TODO-DECISION: Above describes input & output (?) data sets. How relates to s/w modes?
            %   
            % PROPOSAL: Merge LFR and TDS loops.

            inputDatasetLevelList     = {'L1R'};
            inputDashEList            = {'-E'};
            swmSuffixList             = {''};
            swmPurposeAmendmList      = {''};
            if SETTINGS.get_fv('SW_MODES.L1_LFR_TDS_ENABLED')
                inputDatasetLevelList{end+1} = 'L1';
                inputDashEList{end+1}        = '';
                swmSuffixList{end+1}         = '_L1';
                swmPurposeAmendmList{end+1}  = ' EXPERIMENTAL.';
            end



            % Define function which interprets (replaces) specific substrings.            
            % "strmod" = string modify, "g"=global
            strmodg = @(s, iInputLevel) bicas.utils.strrepmany(s, ...
                '<InLvl>',              inputDatasetLevelList{iInputLevel}, ...
                '<I-E>',                inputDashEList{iInputLevel}, ...
                '<SWM suffix>',         swmSuffixList{iInputLevel}, ...
                '<SWM purpose amendm>', swmPurposeAmendmList{iInputLevel});
            
            % Input def that is reused multiple times.
            HK_INPUT_DEF = bicas.swmode_defs.def_input_dataset(...
                'in_hk', 'SOLO_HK_RPW-BIA', 'HK_cdf');

            
            
            LFR_SW_MODE_DATA = struct(...
                'SBMx_SURV',       {'SBM1', 'SBM2', 'SURV', 'SURV'}, ...
                'CWF_SWF',         {'CWF',  'CWF',  'CWF',  'SWF'}, ...
                'modeStr',         {...
                'selective burst mode 1', ...
                'selective burst mode 2', ...
                'survey mode', ...
                'survey mode' ...
                }, ...
                'outputSkeletonVersion', {'09', '09', '09', '09'});
            TDS_SW_MODE_DATA = struct(...
                'CWF_RSWF',              {'CWF', 'RSWF'}, ...
                'outputSkeletonVersion', {'09',  '09'});
            
            

            CUR_INPUT_DEF = bicas.swmode_defs.def_input_dataset('in_cur', 'SOLO_L1_RPW-BIA-CURRENT', 'CUR_cdf');
            
            SwModeList = EJ_library.utils.empty_struct([0,1], ...
                'prodFunc', 'cliOption', 'swdPurpose', 'inputsList', 'outputsList');
            for iInputLevel = 1:numel(inputDatasetLevelList)
                
                %==============================================
                % Iterate over the "fundamental" LFR S/W modes
                %==============================================
                for iSwm = 1:length(LFR_SW_MODE_DATA)
                    strmod = @(s) strrep(strrep(strrep(strmodg(s, iInputLevel), ...
                        '<SBMx/SURV>',  LFR_SW_MODE_DATA(iSwm).SBMx_SURV), ...
                        '<C/SWF>',      LFR_SW_MODE_DATA(iSwm).CWF_SWF), ...
                        '<mode str>',   LFR_SW_MODE_DATA(iSwm).modeStr);

                    SCI_INPUT_DEF = bicas.swmode_defs.def_input_dataset(...
                        'in_sci', ...
                        strmod('SOLO_<InLvl>_RPW-LFR-<SBMx/SURV>-<C/SWF><I-E>'), ...
                        'SCI_cdf');

                    SCI_OUTPUT_DEF = bicas.swmode_defs.def_output_dataset(...
                        strmod('SOLO_L2_RPW-LFR-<SBMx/SURV>-<C/SWF>-E'), ...
                        strmod('LFR L2 <C/SWF> science electric <mode str> data'), ...
                        strmod('RPW LFR L2 <C/SWF> science electric (potential difference) data in <mode str>, time-tagged'), ...
                        LFR_SW_MODE_DATA(iSwm).outputSkeletonVersion);

                    SwModeList(end+1) = bicas.swmode_defs.def_swmode(...
                        @(InputDatasetsMap, Cal) bicas.proc.produce_L2_LFR(...
                            InputDatasetsMap, ...
                            Cal, ...
                            SCI_INPUT_DEF.datasetId, ...
                            SCI_OUTPUT_DEF.datasetId, ...
                            SETTINGS, L), ...
                        strmod('LFR-<SBMx/SURV>-<C/SWF>-E<SWM suffix>'), ...
                        strmod('Generate <SBMx/SURV> <C/SWF> electric field L2 data (potential difference) from LFR <InLvl> data.<SWM purpose amendm>'), ...
                        [SCI_INPUT_DEF, CUR_INPUT_DEF, HK_INPUT_DEF], [SCI_OUTPUT_DEF]);
                end

                %==============================================
                % Iterate over the "fundamental" TDS S/W modes
                %==============================================
                for iSwm = 1:numel(TDS_SW_MODE_DATA)
                    strmod = @(s) strrep(strmodg(s, iInputLevel), ...
                        '<C/RSWF>', TDS_SW_MODE_DATA(iSwm).CWF_RSWF);

                    SCI_INPUT_DEF = bicas.swmode_defs.def_input_dataset(...
                        'in_sci', ...
                        strmod('SOLO_<InLvl>_RPW-TDS-LFM-<C/RSWF><I-E>'), ...
                        'SCI_cdf');
                    
                    SCI_OUTPUT_DEF = bicas.swmode_defs.def_output_dataset(...
                        strmod('SOLO_L2_RPW-TDS-LFM-<C/RSWF>-E'), ...
                        strmod('LFR L2 <C/RSWF> science electric LF mode data'), ...
                        strmod('RPW TDS L2 <C/RSWF> science electric (potential difference) data in LF mode, time-tagged'), ...
                        TDS_SW_MODE_DATA(iSwm).outputSkeletonVersion);

                    SwModeList(end+1) = bicas.swmode_defs.def_swmode(...
                        @(InputDatasetsMap, Cal) bicas.proc.produce_L2_TDS(...
                            InputDatasetsMap, ...
                            Cal, ...
                            SCI_INPUT_DEF.datasetId, ...
                            SCI_OUTPUT_DEF.datasetId, ...
                            SETTINGS, L), ...
                        strmod('TDS-LFM-<C/RSWF>-E<SWM suffix>'), ...
                        strmod('Generate <C/RSWF> electric field L2 data (potential difference) from TDS LF mode <InLvl> data.<SWM purpose amendm>'), ...
                        [SCI_INPUT_DEF, CUR_INPUT_DEF, HK_INPUT_DEF], [SCI_OUTPUT_DEF]);
                end
            end    % for iInputLevel = 1:numel(inputDatasetLevelList)
            
            
            
            obj.List = SwModeList;
            
            EJ_library.assert.castring_set({obj.List(:).cliOption})
        end    % Constructor



        function swModeInfo = get_sw_mode_info(obj, swModeCliOption)
            i = find(strcmp(swModeCliOption, {obj.List(:).cliOption}));
            EJ_library.assert.scalar(i)
            swModeInfo = obj.List(i);
        end
        
        
        
    end    % methods(Access=public)


    methods(Static, Access=private)
        
        
        
        % NOTE: Name dangerously similar to "bicas.swmode_defs".
        function Def = def_swmode(prodFunc, cliOption, swdPurpose, inputsList, outputsList)
            Def.prodFunc    = prodFunc;
            % NOTE: s/w mode CLI _ARGUMENT_ is not intended to be prefixed by e.g. "--". Variable therefore NOT named *Body.
            Def.cliOption   = cliOption;   
            Def.swdPurpose  = swdPurpose;
            Def.inputsList  = inputsList;
            Def.outputsList = outputsList;
            
            
            
            % ASSERTIONS
            EJ_library.assert.castring_set( {...
                Def.inputsList( :).cliOptionHeaderBody, ...
                Def.outputsList(:).cliOptionHeaderBody })   % Important. Check uniqueness of SIP options.
            EJ_library.assert.castring_set( {...
                Def.inputsList(:).prodFuncInputKey })   % Maybe not really necessary.
            EJ_library.assert.castring_set( {...
                Def.outputsList(:).prodFuncOutputKey })   % Maybe not really necessary.
            
            bicas.swmode_defs.assert_SW_mode_CLI_option(Def.cliOption)
            bicas.swmode_defs.assert_text(              Def.swdPurpose)
            
            EJ_library.assert.castring_set({...
                Def.inputsList(:).cliOptionHeaderBody, ...
                Def.outputsList(:).cliOptionHeaderBody})
        end

        
        
        function Def = def_input_dataset(cliOptionHeaderBody, datasetId, prodFuncInputKey)
            % NOTE: No dataset/skeleton version.
            Def.cliOptionHeaderBody = cliOptionHeaderBody;
            Def.prodFuncInputKey    = prodFuncInputKey;
            Def.datasetId           = datasetId;
            
            bicas.swmode_defs.assert_SIP_CLI_option(Def.cliOptionHeaderBody)
            % NOTE: Using the INTERNAL assertion function, not the global one.
            bicas.swmode_defs.assert_DATASET_ID(    Def.datasetId)
        end

        
        
        function Def = def_output_dataset(datasetId, swdName, swdDescription, skeletonVersion)
            Def.cliOptionHeaderBody = 'out_sci';
            Def.prodFuncOutputKey   = 'SCI_cdf';
            Def.swdName             = swdName;
            Def.swdDescription      = swdDescription;
            Def.datasetId           = datasetId;
            Def.datasetLevel        = 'L2';
            Def.skeletonVersion     = skeletonVersion;
            
            bicas.swmode_defs.assert_SW_mode_CLI_option(Def.cliOptionHeaderBody)
            bicas.swmode_defs.assert_text(              Def.swdName)
            bicas.swmode_defs.assert_text(              Def.swdDescription)
            bicas.swmode_defs.assert_DATASET_ID(        Def.datasetId)
            bicas.assert_dataset_level(                 Def.datasetLevel)
            bicas.assert_skeleton_version(              Def.skeletonVersion)
        end



        % NOTE: Wrapper around global counterpart.
        function assert_DATASET_ID(datasetId)
            % PROPOSAL: Use classification function for DATASET_ID instead.
            
            bicas.assert_DATASET_ID(datasetId)
            
            % ASSERTION: Only using SOLO_* DATASET_IDs.
            assert(strcmp('SOLO_', datasetId(1:5)))
        end
        
        

        % Assert that string contains human-readable text.
        function assert_text(str)
            EJ_library.assert.castring_regexp(str, '.* .*')
            EJ_library.assert.castring_regexp(str, '[^<>]*')
        end
        
        
        
        function assert_SW_mode_CLI_option(swModeCliOption)
            EJ_library.assert.castring_regexp(swModeCliOption, bicas.swmode_defs.SW_MODE_CLI_PARAMETER_REGEX)
        end
        
        

        % NOTE: Really refers to "option body".
        function assert_SIP_CLI_option(sipCliOptionBody)
            EJ_library.assert.castring_regexp(sipCliOptionBody, bicas.swmode_defs.SIP_CLI_OPTION_BODY_REGEX)
        end
        
        
        
    end    % methods(Static, Access=private)

end
