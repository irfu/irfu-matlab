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
% Reasons for using a singleton class (instead of static methods:
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
%
% PROPOSAL: More validation.
%   PROPOSAL: Check that data types are unique.
%       NOTE: Requires access to the lists.
%
% PROPOSAL: Use (nested) function to set every input in produce_inputs_constants.
%     PROPOSAL: Same for produce_inputs_constants.
%     PROPOSAL: Use struct statement instead.
%        CON: Does not make use of the similarities between assignments.
%        CON: Want to "extract values from a table".
%     PROPOSAL: Use functions to produce equivalent S/W modes for different input dataset versions (V01-->V02, V02-->V02).
%
% PROPOSAL: Create records (structs) via internal function instead of manually specifying fields.
%   Ex: produce_sw_modes_constants, produce_inputs_constants (partially implemented), produce_outputs_constants.
%   PRO: Forces identical structs: same field names, same set of fields.
%   PRO: Easier to change field names.
%   PRO: Can include assertions(?) on assignment.
%   PRO: More compact assignment code.
%   PRO: Easier to create struct array (where now uses cell array).
%
% PROPOSAL: Use arrays of structs instead of cells.
%    PRO: Forces the use of the same struct fields.
%    NOTE: Would need to create new version of "select_cell_array_structs" that works on arrays instead.
%
% PROPOSAL: Change name to ~dm_constants.
%   NOTE: Should then get rid of BICAS_ROOT_PATH first.
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
            INITIAL_RELEASE_DATE_STR = '2017-02-22';
            INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            
            [obj.INPUTS_INFO_LIST,  obj.INPUTS_PDIDS_LIST]  = bicas.constants.produce_inputs_constants();
            [obj.OUTPUTS_INFO_LIST, obj.OUTPUTS_PDIDS_LIST] = bicas.constants.produce_outputs_constants(INITIAL_RELEASE_DATE_STR, INITIAL_RELEASE_MODIFICATION_STR);          
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
            
            % The RCS ICD, iss2rev2, section 5.3 seems (ambiguous) to imply this regex for CLI S/W mode parameters.
            SW_MODE_CLI_PARAMETER_REGEX = '^[A-Za-z][\w-]+$';   % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.

            % The RCS ICD, iss2rev2, section 3.2.3 only permits these characters (and only lowercase).
            INPUT_CLI_PARAMETER_NAME_PERMITTED_CHARACTERS = 'abcdefghijklmnopqrstuvxyz0123456789_';
            
            %==========================
            % Iterate over input types
            %==========================
            for iInput = 1:length(obj.INPUTS_INFO_LIST)
                cliParameter = obj.INPUTS_INFO_LIST{iInput}.CLI_PARAMETER;
                
                % NOTE: Implicitly checks that cliParameter does NOT begin with "--".
                disallowedCharsFound = setdiff(cliParameter, INPUT_CLI_PARAMETER_NAME_PERMITTED_CHARACTERS);
                if ~isempty(disallowedCharsFound)
                    error('BICAS:constants:Assertion:IllegalConfiguration', ...
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
                    error('BICAS:constants:Assertion:IllegalConfiguration', ...
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
        %    .ID               : S/W mode ID. Used to identify the mode internally (in particular for hardcoded constants
        %                        in data_manager).
        %                        Has about the same purpose as CLI_PARAMETER but is separate so that CLI_PARAMETER
        %                        values/constants can be easily modified, whereas ID values are tied to hardcoded
        %                        constants in data_manager which are harder to modify.
        %    .OUTPUT_PDID_LIST : A cell array of PDIDs. Effectively an array of pointers to (1) the output constants, and (2)
        %                        indirectly to the input constants through data_manager.get_elementary_input_PDIDs.
        function swModesInfoList = produce_sw_modes_constants()
            
            swModesInfoList = {};
            
            %=====
            % LFR 
            %=====
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SBM1-CWF-E_V01-V02';
            SwModeInfo.ID               = 'LFR-SBM1-CWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SBM1-CWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SBM1-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SBM1-CWF-E_V02-V02';
            SwModeInfo.ID               = 'LFR-SBM1-CWF-E_V02-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SBM1-CWF_V02', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SBM1-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SBM2-CWF-E_V01-V02';
            SwModeInfo.ID               = 'LFR-SBM2-CWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SBM2-CWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SBM2-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SBM2-CWF-E_V02-V02';
            SwModeInfo.ID               = 'LFR-SBM2-CWF-E_V02-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SBM2-CWF_V02', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SBM2-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SURV-CWF-E_V01-V02';
            SwModeInfo.ID               = 'LFR-SURV-CWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SURV-CWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SURV-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SURV-CWF-E_V02-V02';
            SwModeInfo.ID               = 'LFR-SURV-CWF-E_V02-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric field data (potential difference) from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SURV-CWF_V02', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SURV-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SURV-SWF-E_V01-V02';
            SwModeInfo.ID               = 'LFR-SURV-SWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate SWF electric (potential difference) data from LFR';            
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SURV-SWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SURV-SWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'LFR-SURV-SWF-E_V02-V02';
            SwModeInfo.ID               = 'LFR-SURV-SWF-E_V02-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate SWF electric (potential difference) data from LFR';
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_LFR-SURV-SWF_V02', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_LFR-SURV-SWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
            %=====
            % TDS
            %=====
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'TDS-LFM-CWF-E_V01-V02';
            SwModeInfo.ID               = 'TDS-LFM-CWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate CWF electric (potential difference) data from TDS-LFM-CWF';
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_TDS-LFM-CWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_TDS-LFM-CWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'TDS-LFM-RSWF-E_V01-V02';
            SwModeInfo.ID               = 'TDS-LFM-RSWF-E_V01-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V01';
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_TDS-LFM-RSWF_V01', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_TDS-LFM-RSWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            SwModeInfo = [];
            SwModeInfo.CLI_PARAMETER    = 'TDS-LFM-RSWF-E_V02-V02';
            SwModeInfo.ID               = 'TDS-LFM-RSWF-E_V02-V02';
            SwModeInfo.SWD_PURPOSE      = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V02';
            SwModeInfo.INPUT_PDID_LIST  = {'L2R_TDS-LFM-RSWF_V02', 'HK_BIA_V02'};
            SwModeInfo.OUTPUT_PDID_LIST = {'L2S_TDS-LFM-RSWF-E_V02'};
            swModesInfoList{end+1} = SwModeInfo;
            
        end



        % Produce constants for all possible INPUT datasets.
        % (independent of how they are associated with S/W modes).
        %
        function [inputsInfoList, einPdidList] = produce_inputs_constants
            % PROPOSAL: Put derivation of .PDID in nested init function.
            
            % NOTE: NESTED function
            function InputInfo = init_input_C_sci(datasetId, skeletonVersionStr)
                InputInfo.CLI_PARAMETER        = 'input_sci';
                InputInfo.DATASET_ID           = datasetId;
                InputInfo.SKELETON_VERSION_STR = skeletonVersionStr;
            end

            inputsInfoList = {};
            
            % NOTE: CLI_PARAMETER = CLI parameter MINUS flag prefix ("--"), i.e. e.g. "sci" instead of "--sci".
            
            %=========
            % BIAS HK
            %=========
            inputsInfoList{end+1} = [];
            inputsInfoList{end}.CLI_PARAMETER        = 'input_hk';
            inputsInfoList{end}.DATASET_ID           = 'ROC-SGSE_HK_RPW-BIA';
            inputsInfoList{end}.SKELETON_VERSION_STR = '01';
            
            inputsInfoList{end+1} = [];
            inputsInfoList{end}.CLI_PARAMETER        = 'input_hk';
            inputsInfoList{end}.DATASET_ID           = 'ROC-SGSE_HK_RPW-BIA';
            inputsInfoList{end}.SKELETON_VERSION_STR = '02';
            
            %=========
            % LFR SCI
            %=========            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '02');   % 1   sample/record
            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '02');   % 1   sample/record
            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '02');   % 1   sample/record
            
            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '01');   % 1 snapshot/record
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '02');   % 1 snapshot/record(!). Adds zVar SAMP_DTIME
            
            %=========
            % TDS SCI
            %=========
            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-CWF',  '01');  % 1   sample/record
            
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '01');  % 1   sample/record.  Adds two zVariables: SAMP_DTIME, SAMPS_PER_CH
            inputsInfoList{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '02');  % 1 snapshot/record

            

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
        % TODO: Set SWD_LEVEL automatically?!
        % PROPOSAL: initialRelaseDateStr, initialRelaseModificationStr as "file-global" constants.
            
            CLI_PARAMETER_SCI_NAME = 'output_sci';
            
            outputsInfoList = {};
            
            % -------- LFR --------
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s CWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 1, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s CWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s CWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s SWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
            
            % -------- TDS --------

            outputsInfoList{end+1} = [];            
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'TDS L2s CWF science electric data in low frequency mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'TDS L2s RSWF science electric data in low frequency mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = initialRelaseDateStr;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = initialRelaseModificationStr;


            % Add one field ".PDID" to every struct above!
            % 
            % Put together PDIDs (used in data_manager).
            % See data_manager for definition.
            eoutPdidList = {};
            for i = 1:length(outputsInfoList)
                outputsInfoList{i}.PDID = bicas.constants.construct_PDID(outputsInfoList{i}.DATASET_ID, outputsInfoList{i}.SKELETON_VERSION_STR);
                eoutPdidList{i} = outputsInfoList{i}.PDID;
            end
        end
        
        
        
        % Construct a PDID derived from a dataset ID and skeleton version (a string shorter than the similar
        % corresponding official strings, e.g. ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E_V01 --> L2S-TDS-LFM-RSWF-E_V01).
        %
        function pdid = construct_PDID(datasetId, skeletonVersionStr)
        
            %pdid = [datasetId, '_V', skeletonVersionStr];
            
            datasetIdShortened = regexprep(datasetId,          '^ROC-SGSE_', '',  'once');
            datasetIdShortened = regexprep(datasetIdShortened, '_RPW-',      '_', 'once');                
            pdid = [datasetIdShortened, '_V', skeletonVersionStr];
        end

    end % methods(Static, Access=private)
    
end   % classdef
