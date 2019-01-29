% Constants - Singleton class for global constants used by BICAS.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software. Set up as a ~singleton handle class.
% Also contains validation code and functions for more convenient access.
%
%
%
% IMPLEMENTATION NOTE
% ===================
% Reasons for using a singleton class (instead of static methods and static class variables):
% 1) Can use properties/instance variables for "caching" values. Do not want to use persistent variables since they
% cause trouble when testing. NOTE: There are no proper static variables in MATLAB.
% 2) Can split up (structure, organize) configuration and validation code in methods.
% 3) The constructor can be used as initialization code which must be run before using the class/constants.
%
%
%
% ~"BUG"?/NOTE: The current implementation contains a minor error of thought(?): It contains an array with data for every
% possible output format. There, every output format is associated with "release data" (required for the S/W descriptor).
% This "release data" should possibly(?) be associated with every S/W mode instead.
%
classdef constants < handle
%
% PROPOSAL: Include SW root path?! How?
%   PRO: Needs to be universally accessible.
%   CON: Not hardcoded. ==> Mixes code with manually set/hardcoded constants.
%   QUESTION: Is there anything analogous? output dir?
%   PROPOSAL: Some functionality for setting "properties", in reality named ~global variables as key-value pairs. cf OVT.
%   PROPOSAL: Handle "manually" through function parameters.
%   TODO-NEED-INFO: There is a ROC-defined environment variable for this?
%
% PROPOSAL: Get rid of BICAS_ROOT_PATH somehow. Does not fit in.
%   PROPOSAL: Move to SETTINGS?
%   PROPOSAL: Use as regular variable submitted through arguments.
%
% PROPOSAL: More validation.
%   PROPOSAL: Check that data types are unique.
%       NOTE: Requires access to the lists.
%
% PROPOSAL: Use functions to produce equivalent S/W modes for different input dataset versions (V01-->V02, V02-->V02).
%
% PROPOSAL: Use arrays of structs instead of cells.
%    PRO: Forces the use of the same struct fields.
%    NOTE: Would need to create new version of "select_cell_array_structs" that works on arrays instead.
%
% PROPOSAL: Change name to which is more specific. Current name is too generic.
%   PROPOSAL: ~dm_constants.
%       NOTE: Should then get rid of BICAS_ROOT_PATH first.
%   PROPOSAL: ~datasets_modes_constants
% PROPOSAL: Split into code for
%   (1) storing the data, and
%   (2) initializing the data
%
% PROPOSAL: Use containers.Map
%   NOTE: Only suitable for structs with exactly one field which is always unique and is "always" used for look-up.
%   PROPOSAL: S/w modes list.
%###################################################################################################################

    properties(Access=public)
        BICAS_ROOT_PATH
        
        SW_MODES_INFO_LIST   % Information associated with S/W modes.
        
        INPUTS_INFO_LIST     % Information associated with input  datasets.
        OUTPUTS_INFO_LIST    % Information associated with output datasets.
        
        INPUTS_PDIDS_LIST
        OUTPUTS_PDIDS_LIST
    end

    properties(Access=private)
        ALL_DATASET_IDS_LIST    % Collect alla known dataset IDs. Useful for assertions.
    end
    
    %###################################################################################################################
    
    methods(Access=public)
        
        % Constructor
        function obj = constants(bicasRootPath)            
            
            obj.BICAS_ROOT_PATH = bicasRootPath;

            % These two values exist in "settings" in principle, but that is just for as long as there has been no
            % official release. After first release, then the two sets should start diverging.
            INITIAL_RELEASE_DATE_STR         = '2018-01-23';
            INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            
            [obj.INPUTS_INFO_LIST,  obj.INPUTS_PDIDS_LIST]  = bicas.constants.produce_inputs_constants();
            [obj.OUTPUTS_INFO_LIST, obj.OUTPUTS_PDIDS_LIST] = bicas.constants.produce_outputs_constants(...
                INITIAL_RELEASE_DATE_STR, INITIAL_RELEASE_MODIFICATION_STR);          
            obj.SW_MODES_INFO_LIST                          = bicas.constants.produce_sw_modes_constants();
            
            
            
            % Extract list (cell array) of unique dataset IDs for input and output datasets.
            obj.ALL_DATASET_IDS_LIST = unique(cellfun(@(s) ({s.DATASET_ID}), [obj.OUTPUTS_INFO_LIST, obj.INPUTS_INFO_LIST])');
                        
            obj.validate
        end
        
        

        function assert_dataset_ID(obj, datasetId)
        % Assert that argument is a valid dataset ID.
        
            if ~ismember(datasetId, obj.ALL_DATASET_IDS_LIST)
                error('BICAS:constants:Assertion', '"%s" is not a valid dataset ID.', datasetId)
            end
        end
        
        function assert_sw_mode_ID(obj, swModeId)
            
            for iMode = 1:length(obj.SW_MODES_INFO_LIST)
                if strcmp(obj.SW_MODES_INFO_LIST{iMode}.ID, swModeId)
                    return
                end
            end
            error('BICAS:constants:Assertion', '"%s" is not a valid S/W mode ID', swModeId)
        end

        function assert_EIn_PDID(obj, einPdid)
            
            for iInput = 1:length(obj.INPUTS_INFO_LIST)
                if strcmp(obj.INPUTS_INFO_LIST{iInput}.PDID, einPdid)
                    return
                end
            end
            error('BICAS:constants:Assertion', '"%s" is not a valid EIn PDID', einPdid)
        end
%         
%         function assert_EOut_PDID(obj, eoutPdid)
%             
%             for i=1:length(obj.OUTPUTS_INFO_LIST)
%                 if strcmp(obj.OUTPUTS_INFO_LIST{i}.PDID, eoutPdid)
%                     return
%                 end
%             end
%             error('BICAS:constants:Assertion', '"%s" is not a valid EOut PDID', eoutPdid)
%         end
    end   % methods(Access=public)
    
    %###################################################################################################################
    
    methods(Access=private)

        % Any code for double-checking the validity of hardcoded constants.
        function validate(obj)
            
            % The RCS ICD, iss2rev2, section 5.3 seems (ambiguous) to imply this regex for S/W mode CLI parameters.
            SW_MODE_CLI_PARAMETER_REGEX = '^[A-Za-z][\w-]+$';   % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.

            % The RCS ICD, iss2rev2, section 3.2.3 only permits these characters (and only lowercase).
            INPUT_CLI_OPTION_HEADER_SH_PERMITTED_CHARACTERS = 'abcdefghijklmnopqrstuvxyz0123456789_';
            
            %==========================
            % Iterate over input types
            %==========================
            for iInput = 1:length(obj.INPUTS_INFO_LIST)
                cliParameter = obj.INPUTS_INFO_LIST{iInput}.OPTION_HEADER_SH;
                
                % NOTE: Implicitly checks that cliParameter does NOT begin with "--".
                disallowedCharsFound = setdiff(cliParameter, INPUT_CLI_OPTION_HEADER_SH_PERMITTED_CHARACTERS);
                if ~isempty(disallowedCharsFound)
                    error('BICAS:constants:Assertion:IllegalCodeConfiguration', ...
                        'Constants value contains illegal character(s). This indicates a pure configuration bug (hard-coded).');
                end
            end
            
            bicas.utils.assert_strings_unique(obj.INPUTS_PDIDS_LIST)
            bicas.utils.assert_strings_unique(obj.OUTPUTS_PDIDS_LIST)            
            
            swModeCliParameterList = cellfun(@(s) ({s.CLI_PARAMETER}), obj.SW_MODES_INFO_LIST);
            swModeIdList           = cellfun(@(s) ({s.ID           }), obj.SW_MODES_INFO_LIST);
            bicas.utils.assert_strings_unique(swModeCliParameterList);
            bicas.utils.assert_strings_unique(swModeIdList);
            
            % ASSERTION: CONSTANTS.SW_MODES_INFO_LIST{i}.CLI_PARAMETER matches validation regexp.
            for iMode = 1:length(obj.SW_MODES_INFO_LIST)
                cliParameter = obj.SW_MODES_INFO_LIST{iMode}.CLI_PARAMETER;
                
                if isempty(regexp(cliParameter, SW_MODE_CLI_PARAMETER_REGEX, 'once'))
                    error('BICAS:constants:Assertion:IllegalCodeConfiguration', ...
                        'Illegal S/W mode CLI parameter definition. This indicates a pure (hard-coded) configuration bug.');
                end
            end
            
            % NOTE: Check that combinations of dataset_ID and SKELETON_VERSION_STR are unique.
            % Implemented by merging strings and checking for unique strings.
            % Is strictly speaking very slightly unsafe; could get false negatives.
            datasetIdVersionList = cellfun( ...
                @(x) ({[x.DATASET_ID, '_V', x.SKELETON_VERSION_STR]}), ...
                [obj.OUTPUTS_INFO_LIST, obj.INPUTS_INFO_LIST]   );
            bicas.utils.assert_strings_unique(datasetIdVersionList)
            
        end

    end   % methods(Access=private)
    
    %###################################################################################################################
    
    methods(Static, Access=private)

        % Define the S/W modes and their associated metadata. 
        % The S/W modes defined here are the only ones which "officially" exist and the only ones which can be used at
        % any given time. The choices here influence (at least) the required CLI arguments and the S/W descriptor.
        %
        % swModesInfoList : cell array of structs
        %    .CLI_PARAMETER    : Is used as CLI parameter to identify the S/W mode.
        %                        NOTE: This is not necessarily to regard as a "CLI option" as defined in "parse_CLI_options".
        %    .ID               : S/W mode ID. Used to identify the mode internally (in particular for hardcoded constants
        %                        in data_manager).
        %                        Has about the same purpose as CLI_PARAMETER but is separate so that CLI_PARAMETER
        %                        values/constants can be easily modified, whereas ID values are tied to hardcoded
        %                        constants in data_manager which are harder to modify.
        %    .OUTPUT_PDID_LIST : A cell array of PDIDs. Effectively an array of pointers to (1) the output constants, and (2)
        %                        indirectly to the input constants through data_manager.get_elementary_input_PDIDs.
        function swModesInfoList = produce_sw_modes_constants()
            % PROPOSAL: Rename CLI_PARAMETER. CLI_NAME? CLI_ARGUMENT?
            % PROPOSAL: Indent hardcoded s/w modes.
            %
            % ~BUG: SWD_PURPOSE should also reference potential data (not just diffs).
            
            % NOTE: NESTED function
            % SMI = Software Mode Info
            function smi = create_SMI(cliParameter, id, inputPdidList, outputPdidList, swdPurpose)
                smi.CLI_PARAMETER    = cliParameter;
                smi.ID               = id;
                smi.SWD_PURPOSE      = swdPurpose;
                smi.INPUT_PDID_LIST  = inputPdidList;
                smi.OUTPUT_PDID_LIST = outputPdidList;
            end            

            swModesInfoList = {};
            
            %=====
            % LFR 
            %=====
            % SBM1
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM1-CWF-E_V01-V03', ...
            'LFR-SBM1-CWF-E_V01-V03', ...
            {'L2R_LFR-SBM1-CWF_V01', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM1-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');            
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM1-CWF-E_V02-V03', ...
            'LFR-SBM1-CWF-E_V02-V03', ...
            {'L2R_LFR-SBM1-CWF_V02', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM1-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');            
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM1-CWF-E_V04-V03', ...
            'LFR-SBM1-CWF-E_V04-V03', ...
            {'L1R_LFR-SBM1-CWF-E_V04', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM1-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');            

            % SBM2
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM2-CWF-E_V01-V03', ...
            'LFR-SBM2-CWF-E_V01-V03', ...
            {'L2R_LFR-SBM2-CWF_V01', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM2-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');            
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM2-CWF-E_V02-V03', ...
            'LFR-SBM2-CWF-E_V02-V03', ...
            {'L2R_LFR-SBM2-CWF_V02', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM2-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SBM2-CWF-E_V04-V03', ...
            'LFR-SBM2-CWF-E_V04-V03', ...
            {'L1R_LFR-SBM2-CWF-E_V04', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SBM2-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');

            % SURV-CWF
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SURV-CWF-E_V01-V03', ...
            'LFR-SURV-CWF-E_V01-V03', ...
            {'L2R_LFR-SURV-CWF_V01', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SURV-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');
            %swModesInfoList{end+1} = create_SMI(...
            %'LFR-SURV-CWF-E_V02-V03', ...
            %'LFR-SURV-CWF-E_V02-V03', ...
            %{'L2R_LFR-SURV-CWF_V02', 'HK_BIA_V02'}, ...
            %{'L2S_LFR-SURV-CWF-E_V03'}, ...
            %'Generate CWF electric field data (potential difference) from LFR');
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SURV-CWF-E_V04-V03', ...
            'LFR-SURV-CWF-E_V04-V03', ...
            {'L1R_LFR-SURV-CWF-E_V04', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SURV-CWF-E_V03'}, ...
            'Generate CWF electric field data (potential difference) from LFR');
            %swModesInfoList{end+1} = create_SMI(...
            %'SOLO_LFR-SURV-CWF-E_V01-V01', ...
            %'SOLO_LFR-SURV-CWF-E_V01-V01', ...
            %{'SOLO_L1R_LFR-SURV-CWF-E_V01', 'HK_BIA_V02'}, ...
            %{'SOLO_L2_LFR-SURV-CWF-E_V01'}, ...
            %'Generate CWF electric field data (potential difference) from LFR');   % RODP TEST

            % SURV-SWF
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SURV-SWF-E_V01-V03', ...
            'LFR-SURV-SWF-E_V01-V03', ...
            {'L2R_LFR-SURV-SWF_V01', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SURV-SWF-E_V03'}, ...
            'Generate SWF electric (potential difference) data from LFR');
            %swModesInfoList{end+1} = create_SMI(...
            %'LFR-SURV-SWF-E_V02-V02', ...
            %'LFR-SURV-SWF-E_V02-V02', ...
            %{'L2R_LFR-SURV-SWF_V02', 'HK_BIA_V02'}, ...
            %{'L2S_LFR-SURV-SWF-E_V03'}, ...
            %'Generate SWF electric (potential difference) data from LFR');
            swModesInfoList{end+1} = create_SMI(...
            'LFR-SURV-SWF-E_V04-V03', ...
            'LFR-SURV-SWF-E_V04-V03', ...
            {'L1R_LFR-SURV-SWF-E_V04', 'HK_BIA_V02'}, ...
            {'L2S_LFR-SURV-SWF-E_V03'}, ...
            'Generate SWF electric (potential difference) data from LFR');

            %=====
            % TDS
            %=====
            
            if 0
                % Modes disabled since they do not work yet.
                % ==> Keep corresponding processing, input & output datasets.
                % NOTE: Not updated to latest output skeletons.
                
                % CWF
                swModesInfoList{end+1} = create_SMI(...
                    'TDS-LFM-CWF-E_V01-V02', ...
                    'TDS-LFM-CWF-E_V01-V02', ...
                    {'L2R_TDS-LFM-CWF_V01', 'HK_BIA_V02'}, ...
                    {'L2S_TDS-LFM-CWF-E_V02'}, ...
                    'Generate CWF electric (potential difference) data from TDS-LFM-CWF');
                
                % RSWF
                swModesInfoList{end+1} = create_SMI(...
                    'TDS-LFM-RSWF-E_V01-V02', ...
                    'TDS-LFM-RSWF-E_V01-V02', ...
                    {'L2R_TDS-LFM-RSWF_V01', 'HK_BIA_V02'}, ...
                    {'L2S_TDS-LFM-RSWF-E_V02'}, ...
                    'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V01');
                swModesInfoList{end+1} = create_SMI(...
                    'TDS-LFM-RSWF-E_V02-V02', ...
                    'TDS-LFM-RSWF-E_V02-V02', ...
                    {'L2R_TDS-LFM-RSWF_V02', 'HK_BIA_V02'}, ...
                    {'L2S_TDS-LFM-RSWF-E_V02'}, ...
                    'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V02');
            end

        end



        % Produce constants for all possible INPUT datasets.
        % (independent of how they are associated with S/W modes).
        %
        function [inputsInfoList, einPdidList] = produce_inputs_constants
            % PROPOSAL: Put derivation of .PDID in nested init function.
            
            % NOTE: NESTED function
            % II = Input Info
            function InputInfo = create_II(datasetId, skeletonVersionStr)
                % PROPOSAL: skeletonVersionStr first argument.
                %   PRO: Rows line up automatically.
                InputInfo.OPTION_HEADER_SH     = 'input_sci';
                InputInfo.DATASET_ID           = datasetId;
                InputInfo.SKELETON_VERSION_STR = skeletonVersionStr;
            end

            inputsInfoList = {};
            
            % NOTE: OPTION_HEADER_SH = CLI option header MINUS option prefix ("--"), i.e. e.g. "sci" instead of "--sci".
            %       SH = short
            
            %=========
            % BIAS HK
            %=========
            %ii = [];    % II = input info
            %ii.OPTION_HEADER_SH     = 'input_hk';
            %ii.DATASET_ID           = 'ROC-SGSE_HK_RPW-BIA';
            %ii.SKELETON_VERSION_STR = '01';
            %inputsInfoList{end+1} = ii;
            
            ii = [];
            ii.OPTION_HEADER_SH     = 'input_hk';
            ii.DATASET_ID           = 'ROC-SGSE_HK_RPW-BIA';
            ii.SKELETON_VERSION_STR = '02';
            inputsInfoList{end+1} = ii;
            
            % RODP TEST
            %ii = [];
            %ii.OPTION_HEADER_SH     = 'input_hk';
            %ii.DATASET_ID           = 'SOLO_HK_RPW-BIA';
            %ii.SKELETON_VERSION_STR = '01';
            %inputsInfoList{end+1} = ii;
            
            %=========
            % LFR SCI
            %=========
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '02');   % 1   sample/record
            inputsInfoList{end+1} = create_II('ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E', '04'); % 1   sample/record
            
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '02');   % 1   sample/record
            inputsInfoList{end+1} = create_II('ROC-SGSE_L1R_RPW-LFR-SBM2-CWF-E', '04'); % 1   sample/record
            
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '01');   % 1 snapshot/record
            %inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '02');   % 1   sample/record
            inputsInfoList{end+1} = create_II('ROC-SGSE_L1R_RPW-LFR-SURV-CWF-E', '04');   % 1 snapshot/record
            %inputsInfoList{end+1} = create_II(    'SOLO_L1R_RPW-LFR-SURV-CWF-E', '01');   % RODP TEST
            
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '01');   % 1 snapshot/record
            %inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '02');   % 1 snapshot/record(!). Adds zVar SAMP_DTIME
            inputsInfoList{end+1} = create_II('ROC-SGSE_L1R_RPW-LFR-SURV-SWF-E', '04'); % 1 snapshot/record(?).
            
            %=========
            % TDS SCI
            %=========
            
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-TDS-LFM-CWF',  '01');  % 1   sample/record
            
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '01');  % 1   sample/record. Adds two zVariables: SAMP_DTIME, SAMPS_PER_CH
            inputsInfoList{end+1} = create_II('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '02');  % 1 snapshot/record

            

            % Add one field ".PDID" to every struct above!
            % 
            % Put together PDIDs (used in data_manager).
            % See data_manager for definition.
            einPdidList = {};
            for i = 1:length(inputsInfoList)
                inputsInfoList{i}.PDID = bicas.constants.construct_PDID(inputsInfoList{i}.DATASET_ID, inputsInfoList{i}.SKELETON_VERSION_STR);
                einPdidList{i} = inputsInfoList{i}.PDID;
            end
        end



        % Produce constants for all possible OUTPUT datasets
        % (independent of how they are associated with S/W modes).
        %
        % ARGUMENTS
        % =========
        % initialRelaseDateStr, initialRelaseModificationStr : For now, values used for all outputs. Should
        %                                                      ideally(?) be set individually for every output.
        function [outputsInfoList, eoutPdidList] = produce_outputs_constants(initialRelaseDateStr, initialRelaseModificationStr)
            
            % NOTE: NESTED function
            % IO = output info
            function OutputInfo = create_OI(datasetId, skeletonVersionStr, swdName, swdDescription)
                % PROPOSAL: skeletonVersionStr first argument.
                %   PRO: Rows line up automatically.
                OutputInfo.SWD_OUTPUT_FILE_IDENTIFIER = 'output_sci';
                OutputInfo.DATASET_ID                 = datasetId;
                OutputInfo.SKELETON_VERSION_STR       = skeletonVersionStr;
                OutputInfo.SWD_NAME                   = swdName;
                OutputInfo.SWD_DESCRIPTION            = swdDescription;
                OutputInfo.SWD_LEVEL                  = 'L2S';
                OutputInfo.SWD_RELEASE_DATE           = initialRelaseDateStr;
                OutputInfo.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
                
                OutputInfo.PDID = bicas.constants.construct_PDID(datasetId, skeletonVersionStr);
            end
            
            outputsInfoList = {};

            % -------- LFR --------
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E', '03', ...
                'LFR L2s CWF science electric data in survey mode', ...
            'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 1, time-tagged');
           
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E', '03', ...
                'LFR L2s CWF science electric data in survey mode', ...
            'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged');
            
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E', '03', ...
                'LFR L2s CWF science electric data in survey mode', ...
            'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged');
            
            %outputsInfoList{end+1} = create_OI('SOLO_L2_RPW-LFR-SURV-CWF-E', '01', ...
            %    'LFR L2 CWF science electric data in survey mode', ...
            %'RPW LFR L2 CWF science electric (potential difference) data in survey mode, time-tagged');   % RODP TEST
            
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E', '03', ...
                'LFR L2s SWF science electric data in survey mode', ...
            'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged');
            
            % -------- TDS --------
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E', '02', ...
                'TDS L2s CWF science electric data in low frequency mode', ...
            'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged');
            
            outputsInfoList{end+1} = create_OI('ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E', '02', ...
                'TDS L2s RSWF science electric data in low frequency mode', ...
            'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged');


        
            % Compile list of Elementary output PDIDs.
            eoutPdidList = {};
            for i = 1:length(outputsInfoList)
                eoutPdidList{i} = outputsInfoList{i}.PDID;
            end
        end
        
        
        
        % Construct a PDID derived from a dataset ID and skeleton version (a string shorter than the similar
        % corresponding official strings, e.g. ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E_V01 --> L2S-TDS-LFM-RSWF-E_V01).
        %
        % NOTE: Has to work sensibly for both ROC-SGSE and RODP/SOLO dataset IDs.
        function pdid = construct_PDID(datasetId, skeletonVersionStr)
        
            %pdid = [datasetId, '_V', skeletonVersionStr];
            
            datasetIdShortened = regexprep(datasetId,          '^ROC-SGSE_', '',  'once');
            datasetIdShortened = regexprep(datasetIdShortened, '_RPW-',      '_', 'once');                
            pdid = [datasetIdShortened, '_V', skeletonVersionStr];
        end

    end % methods(Static, Access=private)
    
end   % classdef
