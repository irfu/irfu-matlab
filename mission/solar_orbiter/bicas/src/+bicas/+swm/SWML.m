%
% Singleton class that stores (after having "built" it) an unmodifiable data
% structure that represents which and how s/w modes are CURRENTLY VISIBLE to the
% user. What that data structure contains thus depends on
% -- current pipeline: RODP, ROC-SGSE
% -- whether support for L1 input datasets is enabled or not.
%
% Data here
% -- is intended to (1) add metadata to the production functions for
%       (a) the caller interface
%       (b) potentially future help text
%       (c) the s/w descriptor
% -- contains variables to make it possible to match information for input and
%    output datasets here, with that of bicas.proc' production functions.
%
%
% IMPLEMENTATION NOTES
% ====================
% The class essentially consists of one large struct, and a constructor that
% builds it. The large data struct contains many parts which are similar but not
% the same. To do this, much of the data is "generated" with hard-coded strings
% (mostly the same in every iteration), in which specific codes/substrings are
% substituted algorithmically (different in different iterations). To avoid
% mistakes, the code uses a lot of assertions to protect against mistakes, e.g.
% -- algorithmic bugs
% -- mistyped hard-coded info
% -- mistakenly confused arguments with each other.
% Assertions are located at the place where "values are placed in their final
% location".
% 
% NOTE: To implement compatibility with L1 input datasets, the code must be able
% to handle
% -- changing input dataset levels: L1 (unofficial support), L1R (official
%    support). It implements support for L1 input datasets via separate S/W
%    modes.
% 
%
% RATIONALE
% =========
% -- Should decrease the amount of overlapping hard-coded information to e.g.
%    reduce risk of mistakes, reduce manual work when verifying updates.
% -- Having one big, albeit somewhat redundant data structure should make the
%    interface to the rest of BICAS relatively future-proof, in the face of
%    future updates.
% -- Useful for expected future bias current datasets.
% -- Possible need for backward compatibility.
%
%
% NOTE
% ====
% There is no s/w mode for generating VHT datasets. They are therefore not
% represented here.
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
classdef SWML
    % PROPOSAL: Pick SWD name/descriptions from master CDFs.
    % PROPOSAL: Obtain output dataset level from production function metadata?!!
    %
    % PROPOSAL: Always produce all possible s/w modes (both pipelines, incl. L1),
    %           then filter out the undesired ones using internal metadata for
    %           every S/W mode.
    %
    % PROPOSAL: Use PF = prodFunc, production function
    % PROPOSAL: Same input CDF can have multiple DATASET_IDs, but only one is
    %           shown in the s/w descriptor.
    %   PRO: Can handle old datasets with ROG-SGSE DATASET_IDs, and otherwise
    %        only use RODP DATASET_IDs.
    %
    % PROPOSAL: Fieldname change
    %   inputsList  --> inputsArray
    %   outputsList --> outputsArray
    %   NOTE: Likely influences BICAS testing code and pipeline.
    %         Should only be implemented at the right time.
    %
    % TODO-DEC: Which arguments should SWML production functions (function handles in
    %           an instance of bicas.swm.SWML) have?
    %   NOTE: The arguments needed by the underlying production functions
    %         varies, but the arguments returned by bicas.swm.SWML must be the same.
    %   NOTE: produce_L1R_to_L2_LFR/TDS() are used for multiple s/w modes with some
    %         arguments hard-coded differently for different s/w modes (input & output DATASET_IDs).
    %   NOTE: SWM/underlying production functions can receive argument values via
    %       (1) bicas.swm.SWML (constructor), or (2) the call in execute_SWM.
    %   PROPOSAL: All arguments which are known at the time bicas.swm.SWML
    %       constructor is called, should receive values there.
    %       ==> ~As many as possible.
    %       CON: bicas.swm.SWML not really meant to set production function arguments.
    %       CON: Makes bicas.swm.SWML harder to initialize (outside of BICAS).
    %   PROPOSAL: All arguments which are different for different (underlying) production
    %             functions. ==> As few as possible.
    %   Ex: SETTINGS, L, rctDir, NsoTable
    %
    % PROPOSAL: Refactor to separate the hardcoded SWMs from the class.
    %   PRO: Good for testing.
    %   PRO: Good for separating hardcoding/"constants".
    %
    % PROPOSAL: Use classes instead of structs.
    %   PROPOSAL: Use sub-package for classes.
    %   PROPOSAL: Class for s/w mode definition.
    %       PROPOSAL: Name SWM
    %           CON: Too similar to bicas.swm.SWML.
    %   PROPOSAL: Class for input datasets.
    %   PROPOSAL: Class for output datasets.
    %
    % PROPOSAL: Refactor to
    %   Redefine and rename as class SWML = SWM List
    %   Class for single SWM.
    %   Function for hardcoding that initializes SWML.
    %   Create package swm.* including swm.swml, swm.hardcoding, swm.swm (class
    %   for single SWM), etc.
    %   get_SWM_info --> get_SWM
    %   PROPOSAL: Abolish this class. Only need SWM class + object array/list.
    %   PROPOSAL:
    %       bicas.swm.swm
    %       bicas.swm.swml
    %       bicas.swm.InputDataset    ?
    %       bicas.swm.OutputDataset   ?
    %       bicas.swm.utils  ?
    %       bicas.swm.get_SWML (function)
    %
    % PROPOSAL: Better class name
    %   PROPOSAL: SwmSet
    %       PRO: There is no inherent ordering of SWMs.
    %   PROPOSAL: SwmList
    %       CON: List implies ordering.



    % PRIVATE, STATIC, CONSTANTS
%     properties(Constant, GetAccess=private)
%
%         SWM_PURPOSE_AMENDMENT = ' UNOFFICIAL wrt. ROC.';
%     end

    
    
    % PUBLIC, IMMUTABLE
    properties(SetAccess=immutable)
        
        % NOTE: Implicit that it is a list of s/w modes (not in name). Note that
        % it is a public property.
        List
    end
    
    
    
    methods(Access=public)
        
        
        
        % Constructor
        %
        % ARGUMENTS
        % =========
        %
        % IMPLEMENTATION NOTE: The constructor used to be written so that it was
        % easy to disable S/W modes with L2R input datasets (for backward
        % compatibility). That functionality has now been now removed, although
        % the implementation has not been entirely updated to take advantage of
        % this (not simplified of this).
        % 
        function obj = SWML(SETTINGS, L)
            % PROPOSAL: Re-implement (top-level) hard-coded constants by setting
            %   multiple redundant 1D(?) vectors that covers every case. Then
            %   set various cases by assigning constants to many elements using
            %   MATLAB syntax. One index representing: Combination of
            %   DATASET_ID+Skeleton_Version (both pipelines, LFR+TDS, HK+SCI),
            %   every element contains data for that dataset. Must use
            %   combination DATASET_ID+Skeleton_Version to potentially cover old
            %   versions. Manipulate and set multiple elements smoothly by using
            %   vectors for indices.
            %   Ex: Vectors to set: skeletonVersionVector, SBMx_SURV_vector,
            %       CWF_SWF_vector, output dataset level, vectors for
            %       human-readable description string(s) (e.g. modeStr)
            %   Ex: Vectors with indices for dataset in/for: either pipeline,
            %       science or HK, LFR or TDS, latest versions or
            %       backward-compatibility versions.
            %   --
            %   TODO-DEC: Above describes input & output (?) data sets.
            %       How relates to s/w modes?
            %   
            % PROPOSAL: Merge LFR and TDS loops.
            % PROPOSAL: Split constructor into separate functions, one per s/w
            %           mode (or at least the big ones)
            %   L1/L1R --> L2 mode
            %   L2     --> L3

            % Arrays with constants.
            % {1} = S/w modes (science) L1R-->L2
            % {2} = S/w modes (science) L1 -->L2
            inputDatasetLevelList     = {'L1R'};
            inputDashEList            = {'-E'};
            swmSuffixList             = {''};
            swmPurposeAmendmList      = {''};
            if SETTINGS.get_fv('SWM.L1-L2_ENABLED')
                inputDatasetLevelList{end+1} = 'L1';
                inputDashEList{end+1}        = '';
                swmSuffixList{end+1}         = '_L1';
                swmPurposeAmendmList{end+1}  = ' UNOFFICIAL wrt. ROC.';
            end



            % Define function which interprets (replaces) specific substrings.            
            % "strmod" = string modify, "g"=global
            strmodg = @(s, iInputLevel) bicas.utils.strrep_many(s, ...
                '<InLvl>',              inputDatasetLevelList{iInputLevel}, ...
                '<I-E>',                inputDashEList{iInputLevel}, ...
                '<SWM suffix>',         swmSuffixList{iInputLevel}, ...
                '<SWM purpose amendm>', swmPurposeAmendmList{iInputLevel});

            % Input definitions that are reused multiple times.
            HK_INPUT_DEF = bicas.swm.InputDataset(...
                'in_hk', 'SOLO_HK_RPW-BIA', 'HK_cdf');
            CUR_INPUT_DEF = bicas.swm.InputDataset(...
                'in_cur', 'SOLO_L1_RPW-BIA-CURRENT', 'CUR_cdf');
            
            
            
            % Define arrays of data used for generating s/w modes definitions.
            % One component per pair of s/w modes L1/L1R-->L2
            LFR_SWM_DATA = struct(...
                'SBMx_SURV',       {'SBM1', 'SBM2', 'SURV', 'SURV'}, ...
                'CWF_SWF',         {'CWF',  'CWF',  'CWF',  'SWF' }, ...
                'modeStr',         {...
                    'selective burst mode 1', ...
                    'selective burst mode 2', ...
                    'survey mode', ...
                    'survey mode' ...
                }, ...
                'outputSkeletonVersion', {'14', '14', '14', '14'});
            TDS_SWM_DATA = struct(...
                'CWF_RSWF',              {'CWF', 'RSWF'}, ...
                'outputSkeletonVersion', {'14',  '14'});
            
            

            SwmList = bicas.swm.SWM.empty(0, 1);
            % Iterate over [L1R] (one component), or [L1, L1R]...
            for iInputLevel = 1:numel(inputDatasetLevelList)
                
                %==============================================
                % Iterate over the "fundamental" LFR S/W modes
                %==============================================
                for iSwm = 1:length(LFR_SWM_DATA)
                    
                    % Define local string modification function.
                    strmod = @(s) bicas.utils.strrep_many(strmodg(s, iInputLevel), ...
                        '<SBMx/SURV>',  LFR_SWM_DATA(iSwm).SBMx_SURV, ...
                        '<C/SWF>',      LFR_SWM_DATA(iSwm).CWF_SWF, ...
                        '<mode str>',   LFR_SWM_DATA(iSwm).modeStr);

                    SciInputDef = bicas.swm.InputDataset(...
                        'in_sci', ...
                        strmod('SOLO_<InLvl>_RPW-LFR-<SBMx/SURV>-<C/SWF><I-E>'), ...
                        'SCI_cdf');

                    SciOutputDef = bicas.swm.SWML.def_output_dataset(...
                        'out_sci', ...
                        strmod('SOLO_L2_RPW-LFR-<SBMx/SURV>-<C/SWF>-E'), ...
                        'SCI_cdf', ...
                        strmod('LFR L2 <C/SWF> science electric <mode str> data'), ...
                        strmod(...
                            ['RPW LFR L2 <C/SWF> science electric', ...
                            ' (potential difference) data in <mode str>,', ...
                            ' time-tagged']), ...
                        LFR_SWM_DATA(iSwm).outputSkeletonVersion);

                    SwmList(end+1) = bicas.swm.SWM(...
                        @(InputDatasetsMap, rctDir, NsoTable) bicas.proc.pf.produce_L1R_to_L2_LFR(...
                            InputDatasetsMap, ...
                            rctDir, ...
                            NsoTable, ...
                            SciInputDef.datasetId, ...
                            SciOutputDef.datasetId, ...
                            SETTINGS, L), ...
                        strmod('LFR-<SBMx/SURV>-<C/SWF>-E<SWM suffix>'), ...
                        strmod(...
                            ['Generate <SBMx/SURV> <C/SWF> electric field', ...
                            ' L2 data (potential difference)', ...
                            ' from LFR <InLvl> data.<SWM purpose amendm>']), ...
                        [SciInputDef, CUR_INPUT_DEF, HK_INPUT_DEF], ...
                        [SciOutputDef]);
                end

                %==============================================
                % Iterate over the "fundamental" TDS S/W modes
                %==============================================
                for iSwm = 1:numel(TDS_SWM_DATA)
                    
                    % Define local string modification function.
                    strmod = @(s) strrep(strmodg(s, iInputLevel), ...
                        '<C/RSWF>', TDS_SWM_DATA(iSwm).CWF_RSWF);

                    SciInputDef = bicas.swm.InputDataset(...
                        'in_sci', ...
                        strmod('SOLO_<InLvl>_RPW-TDS-LFM-<C/RSWF><I-E>'), ...
                        'SCI_cdf');
                    
                    SciOutputDef = bicas.swm.SWML.def_output_dataset(...
                        'out_sci', ...
                        strmod('SOLO_L2_RPW-TDS-LFM-<C/RSWF>-E'), ...
                        'SCI_cdf', ...
                        strmod('LFR L2 <C/RSWF> science electric LF mode data'), ...
                        strmod(...
                            ['RPW TDS L2 <C/RSWF> science electric (potential', ...
                            ' difference) data in LF mode, time-tagged']), ...
                        TDS_SWM_DATA(iSwm).outputSkeletonVersion);

                    SwmList(end+1) = bicas.swm.SWM(...
                        @(InputDatasetsMap, rctDir, NsoTable) bicas.proc.pf.produce_L1R_to_L2_TDS(...
                            InputDatasetsMap, ...
                            rctDir, ...
                            NsoTable, ...
                            SciInputDef.datasetId, ...
                            SciOutputDef.datasetId, ...
                            SETTINGS, L), ...
                        strmod('TDS-LFM-<C/RSWF>-E<SWM suffix>'), ...
                        strmod(...
                            ['Generate <C/RSWF> electric field L2 data', ...
                            ' (potential difference)', ...
                            ' from TDS LF mode <InLvl> data.<SWM purpose amendm>']), ...
                        [SciInputDef, CUR_INPUT_DEF, HK_INPUT_DEF], ...
                        [SciOutputDef]);
                end
            end    % for iInputLevel = 1:numel(inputDatasetLevelList)
            
            
            
            if SETTINGS.get_fv('SWM.L2-L2_CWF-DSR_ENABLED')
                SciInputDef = bicas.swm.InputDataset(...
                    'in_sci', ...
                    'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
                    'OSR_cdf');
                
                SciOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_dsr', ...
                    'SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND', ...
                    'DSR_cdf', ...
                    'LFR L2 CWF science electric survey mode data, downsampled', ...
                    ['RPW LFR L2 CWF science electric (potential difference)', ...
                    ' data in survey mode, time-tagged, downsampled'], ...
                    '02');
                
                % NOTE: Function handle: Argument rctDir is not used, but is
                % needed for the interface.
                SwmList(end+1) = bicas.swm.SWM(...
                    @(InputDatasetsMap, rctDir, NsoTable) bicas.proc.pf.produce_L2_to_L2_CWF_DSR(...
                        InputDatasetsMap, SETTINGS, L), ...
                    'LFR-SURV-CWF-E-DSR', ...
                    ['Generate downsampled version of LFR L2 SURV CWF', ...
                    ' science electric (potential difference) data.', ...
                    ' NOTE: This is an unofficial s/w mode.'], ...
                    [SciInputDef], ...
                    [SciOutputDef]);
                    
            end
            
            
            
            if SETTINGS.get_fv('SWM.L2-L3_ENABLED')
                SciInputDef = bicas.swm.InputDataset(...
                    'in_sci', ...
                    'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
                    'LFR-SURV-CWF-E_cdf');
                
                EfieldOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_efield', ...
                    'SOLO_L3_RPW-BIA-EFIELD', ...
                    'EFIELD_OSR_cdf', ...
                    'BIAS L3 science electric field vector data', ...
                    'RPW BIAS L3 science electric field vector data, time-tagged', ...
                    '04');
                
                ScpotOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_scpot', ...
                    'SOLO_L3_RPW-BIA-SCPOT', ...
                    'SCPOT_OSR_cdf', ...
                    'BIAS L3 science spacecraft potential data', ...
                    'RPW BIAS L3 science spacecraft potential data, time-tagged', ...
                    '04');

                DensityOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_density', ...
                    'SOLO_L3_RPW-BIA-DENSITY', ...
                    'DENSITY_OSR_cdf', ...
                    'BIAS L3 science plasma density data', ...
                    'RPW BIAS L3 science plasma density data, time-tagged', ...
                    '04');

                EfieldDsrOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_efield_dsr', ...
                    'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
                    'EFIELD_DSR_cdf', ...
                    'BIAS L3 downsampled science electric field vector data', ...
                    ['RPW BIAS L3 downsampled science electric', ...
                    ' field vector data, time-tagged'], ...
                    '04');
                
                ScpotDsrOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_scpot_dsr', ...
                    'SOLO_L3_RPW-BIA-SCPOT-10-SECONDS', ...
                    'SCPOT_DSR_cdf', ...
                    'BIAS L3 downsampled science spacecraft potential data', ...
                    ['RPW BIAS L3 downsampled science spacecraft', ...
                    ' potential data, time-tagged'], ...
                    '04');

                DensityDsrOutputDef = bicas.swm.SWML.def_output_dataset(...
                    'out_density_dsr', ...
                    'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS', ...
                    'DENSITY_DSR_cdf', ...
                    'BIAS L3 downsampled science plasma density data', ...
                    ['RPW BIAS L3 downsampled science plasma', ...
                    ' density data, time-tagged'], ...
                    '04');

                % NOTE: Function handle: Arguments rctDir, NsoTable are not
                % used, but are needed for the interface.
                SwmList(end+1) = bicas.swm.SWM(...
                    @(InputDatasetsMap, rctDir, NsoTable) (bicas.proc.pf.produce_L2_to_L3(...
                        InputDatasetsMap, ...
                        SETTINGS, L)), ...
                    'BIA-EFIELD-SCPOT-DENSITY', ...
                    ['Generate L3 electric field vector, spacecraft', ...
                    ' potential, and density data', ...
                    ' incl. additional downsampled versions.', ...
                    ' NOTE: This is an unofficial s/w mode.'], ...
                    [SciInputDef], ...
                    [EfieldOutputDef, ...
                     EfieldDsrOutputDef, ...
                     ScpotOutputDef, ...
                     ScpotDsrOutputDef, ...
                     DensityOutputDef, ...
                     DensityDsrOutputDef]);
            end



            obj.List = SwmList;
            
            irf.assert.castring_set({obj.List(:).cliOption})
        end    % Constructor



        function swmInfo = get_SWM_info(obj, swmCliOption)
            i = find(strcmp(swmCliOption, {obj.List(:).cliOption}));
            irf.assert.scalar(i)
            swmInfo = obj.List(i);
        end
        
        
        
    end    % methods(Access=public)


    
    methods(Static, Access=public)



        % Assert that string contains human-readable text.
        function assert_text(str)
            irf.assert.castring_regexp(str, '.* .*')
            irf.assert.castring_regexp(str, '[^<>]*')
        end



        function assert_SWM_CLI_option(swmCliOption)
            irf.assert.castring_regexp(...
                swmCliOption, ...
                bicas.constants.SWM_CLI_OPTION_REGEX)
        end



        % NOTE: Really refers to "option body".
        function assert_SIP_CLI_option(sipCliOptionBody)
            irf.assert.castring_regexp(...
                sipCliOptionBody, ...
                bicas.constants.SIP_CLI_OPTION_BODY_REGEX)
        end



        % NOTE: Wrapper around global counterpart.
        function assert_DATASET_ID(datasetId)
            bicas.assert_BICAS_DATASET_ID(datasetId)
            
            % ASSERTION: Only using SOLO_* DATASET_IDs.
            [sourceName, ~, ~] = solo.adm.disassemble_DATASET_ID(datasetId);
            assert(strcmp(sourceName, 'SOLO'))
        end
        


    end    % methods(Static, Access=public)



    methods(Static, Access=private)



        function Def = def_output_dataset(...
                cliOptionHeaderBody, datasetId, prodFuncOutputKey, ...
                swdName, swdDescription, skeletonVersion)
            
            [~, datasetLevel, ~] = solo.adm.disassemble_DATASET_ID(datasetId);
            
            Def.cliOptionHeaderBody = cliOptionHeaderBody;
            Def.datasetId           = datasetId;
            Def.datasetLevel        = datasetLevel;
            
            Def.prodFuncOutputKey   = prodFuncOutputKey;   % Ex: 'SCI_cdf';
            Def.swdName             = swdName;
            Def.swdDescription      = swdDescription;
            Def.skeletonVersion     = skeletonVersion;
            
            bicas.swm.SWML.assert_SWM_CLI_option(Def.cliOptionHeaderBody)
            bicas.swm.SWML.assert_text(              Def.swdName)
            bicas.swm.SWML.assert_text(              Def.swdDescription)
            bicas.swm.SWML.assert_DATASET_ID(        Def.datasetId)
            solo.adm.assert_dataset_level(           Def.datasetLevel)
            bicas.assert_skeleton_version(           Def.skeletonVersion)
        end



    end    % methods(Static, Access=private)

end
