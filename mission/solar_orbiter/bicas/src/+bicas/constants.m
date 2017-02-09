% Constants - Singleton class for global constants used by BICAS.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software. Set up as a ~singleton handle class.
% Also contains validation code and functions for more convenient access.
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% Some constants (1) correspond exactly to fields in the S/W (JSON) descriptor, and
% (2) are unlikely to be used for anything else. These are labeled with a prefix "SWD_". Other
% variables which do not have the prefix may also be used for the S/W descriptor too but they are
% probably more unambiguous in their meaning.
%
%
% IMPLEMENTATION NOTE 1
% =====================
% BICAS contains other code which builds a structure corresponding to the S/W descriptor (defined
% in the RCS ICD) from the constants structure here.
% Reasons for NOT putting the S/W descriptor structure inside the constants structure:
% (1) Some of the S/W descriptor variables have vague or misleading names ("name", "dataset versions", "dataset IDs")
%     which would (reasoably) have to be represented by MATLAB variables with the same names.
% (2) Some of the S/W descriptor variables are grouped in a way which
%     does not fit the rest of the code (modes[].outputs.output_XX.release in the S/W descriptor structure).
% (3) Some of the S/W descriptor values are really structure field NAMES, but would be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (4) The constants structure would become dependent on the format of the S/W descriptor structure.
%     The latter might change in the future, or be misunderstood in the present (i.e. be changed).
% (5) Some S/W descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (6) It is easier to add automatic checks on the S/W descriptor in the code that derives it.
%
%
% IMPLEMENTATION NOTE 2
% =====================
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
% PROPOSAL: Redefine file into something that makes "decisions"?
%   PROPOSAL: Method for finding the path to the master CDF for a given output (instead of separate function).
%
% PROPOSAL: Include SW root path?! How?
%    PRO: Needs to be universally accessible.
%    CON: Not hardcoded. ==> Mixes code with manually set constants.
%    CON: No good MATLAB implementation of static variables.
%    QUESTION: Is there anything analogous? output dir?
%    PROPOSAL: Some functionality for setting "properties", in reality named ~global variables as key-value pairs. cf OVT.
%    PROPOSAL: Handle "manually" through function parameters.
%
% PROPOSAL: Add S/W descriptor to this structure by calling get_sw_descriptor?
%    PRO: Can force validation to take place in get_sw_descriptor.
%
% PROPOSAL: More validation.
%   PROPOSAL: Check that master CDFs exist, that paths exist.
%       CON: Makes sense to make that kind of check here?!
%   PROPOSAL: Check that data types are unique.
%       NOTE: Requires access to the lists.
%
% PROPOSAL: Have string identifiers for input/output types identical to what "data_manager" uses?
%    PROPOSAL: Use containers.Map instead of cell arrays.
%       CON: Using dataset IDs as keys creates another double use of that value.
%       CON?: Only useful if it obvious which string to use for look-up. Already have "select_structs" for selecting.
%       PRO: Makes it more convenient to internally have multiple input/output dataset versions.
% PROPOSAL: Use arrays of structs instead of cells.
%    PRO: Forces the use of the same struct fields.
%    NOTE: Would need to create new version of "select_structs" that works on arrays instead.
%
% PROPOSAL: Move get_sw_descriptor function to this class.
%   PRO: It has to do with constants.
%       CON: The function really represents code, not constants.
%   PRO: More natural to incorporate validation.
%   PRO: The result could/should be "cached" here.
% PROPOSAL: Store derived S/W descriptor in constants, even if get_sw_descriptor is not there.
%   CON/PROBLEM: get_sw_descriptor is a function of CONSTANTS which must hence first be initialized. 
%       PROPOSAL: Have get_sw_descriptor be a function of subsets of constants.
%
% PROPOSAL: Change constants to not have any public variables. All constants are accessed through functions which cache
% against private instance variables.
%    PRO: Useful for non-trivial instantiation.
%       Ex: Using external get_sw_descriptor to derive SWD variable INSIDE constants.
%           CON: External function would use constants before it has been properly initialized.
%
% PROPOSAL: Use (nested) function to set every input in produce_inputs_constants. Reduce to one-liners.
%     PROPOSAL: Same for produce_inputs_constants.
%     PROPOSAL: Use struct statement instead.
%        CON: Does not make use of the similarities between assignments.
%        CON: Want to "extract values from a table".
%     PROPOSAL: Use functions to produce equivalent S/W modes for different input dataset versions (V01-->V02, V02-->V02).
%
% PROPOSAL: Convert constant cell arrays of structs to arrays of structs: S/W modes, CDF/PDID inputs/outputs.
%   PRO: Simplifies code that constructs cell arrays of the same struct field in multiple cell structs.
%###################################################################################################################

    properties(Access=public)
        C                    % For miscellaneous minor constants which still might require code to be initialized.
        
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
            
            %-------------------------------------------------------------------------------------
            % Common values
            % -------------
            % Only used INDIRECTLY and only INTERNALLY to set the values of the "real" constants.
            %-------------------------------------------------------------------------------------
            D = [];
            D.INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            D.INITIAL_RELEASE_DATE = '2016-11-17';
            %D.SWD_OUTPUT_RELEASE_VERSION = '01';  % For the S/W descriptor output CDFs' release version. Unknown what a sensible value is.
            
            
            
            C = [];

            C.AUTHOR_NAME              = 'Erik P G Johansson';
            C.AUTHOR_EMAIL             = 'erik.johansson@irfu.se';
            C.INSTITUTE                = 'IRF-U';
            C.BICAS_ROOT_PATH          = bicasRootPath;
            C.MASTER_CDFS_RELATIVE_DIR = 'data';    % Location of master CDF files. Relative to the software directory structure root.
            
            % Value that shows up in EOut dataset GlobalAttributes.Calibration_version.
            % String value.
            C.CALIBRATION_VERSION = '0.1; Only proportionality constants i.e. no voltage offset tables, no transfer functions; No bias currents';
        
            
            
            %==========================================================================================================
            % Various S/W descriptor release data for the entire software (not specific outputs)
            % ----------------------------------------------------------------------------------
            % EXCEPTION TO VARIABLE NAMING CONVENTION: Field names are used for constructing the JSON object struct and
            % can therefore NOT follow variable naming conventions without modifying other code.
            %==========================================================================================================
            C.SWD_IDENTIFICATION.project     = 'ROC-SGSE';
            C.SWD_IDENTIFICATION.name        = 'BICAS';
            C.SWD_IDENTIFICATION.identifier  = 'ROC-SGSE-BICAS';
            C.SWD_IDENTIFICATION.description = 'BIAS Calibration Software (BICAS) which derives the BIAS L2S input signals (plus some) from the BIAS L2R output signals.';
            %
            C.SWD_RELEASE.version      = '0.1.0';
            C.SWD_RELEASE.date         = D.INITIAL_RELEASE_DATE;
            C.SWD_RELEASE.author       = C.AUTHOR_NAME;
            C.SWD_RELEASE.contact      = C.AUTHOR_EMAIL;
            C.SWD_RELEASE.institute    = C.INSTITUTE;
            C.SWD_RELEASE.modification = D.INITIAL_RELEASE_MODIFICATION_STR;
            %
            C.SWD_ENVIRONMENT.executable = 'roc/bicas';
            
            
            
            % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
            C.STDOUT_PREFIX = 'STDOUT: ';
        
            % Parameters influencing how JSON objects are printed with function JSON_object_str.
            C.JSON_OBJECT_STR = [];
            C.JSON_OBJECT_STR.INDENT_SIZE    =  4;
            C.JSON_OBJECT_STR.VALUE_POSITION = 15;
        
            % The epoch for ACQUISITION_TIME.
            % The time in UTC at which ACQUISITION_TIME is [0,0].
            % Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
            % PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
            C.ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];

            C.INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID       = 0;   % Require input CDF Global Attribute "DATASET_ID"       to match the expected value.
            C.INPUT_CDF_ASSERTIONS.STRICT_SKELETON_VERSION = 1;   % Require input CDF Global Attribute "Skeleton_version" to match the expected value.
            C.INPUT_CDF_ASSERTIONS.MATCHING_TEST_ID        = 0;   % Require Test_id to be identical for all input CDF datasets.
            C.OUTPUT_CDF.SET_TEST_ID = 1;            % Set CDF GlobalAttribute "Test_id". ROC DFMD says that it should really be set by ROC.            
            C.OUTPUT_CDF.DATA_VERSION = '01';        % Set CDF GlobalAttribute "Data_version". ROC DFMD says it should be updated in a way which can not be automatized?!!! Set here for now.
            
            C.PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION = 1;
            
            % zVariables which are still empty after copying data into the master CDF assigned a correctly sized array
            % with fill values. This should only be necessary for S/W modes with incomplete processing.
            C.OUTPUT_CDF.EMPTY_ZVARIABLES_SET_TO_FILL = 0;
            
            C.LOGGING.MAX_UNIQUES_PRINTED = 5;    % When logging contents of matrix/vector, maximum number of unique values printed before switching to shorter representation (min-max range)
            C.LOGGING.IRF_LOG_LEVEL = 'notice';   % Log level for "irf.log".
            
            
            
            %=====================================================================
            % Define constants relating to interpreting LFR datasets
            % ------------------------------------------------------
            % F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz. Names are LFR's naming.
            %=====================================================================
            C.LFR = [];
            C.LFR.F0 = 24576;  % = 6 * 4096
            C.LFR.F1 =  4096;
            C.LFR.F2 =   256;
            C.LFR.F3 =    16;
        
            %========================================================
            % Constants for how the "simple demuxer" calibrates data
            %========================================================
            C.SIMPLE_DEMUXER = [];
            C.SIMPLE_DEMUXER.ALPHA           = 1/17;
            C.SIMPLE_DEMUXER.BETA            =    1;
            C.SIMPLE_DEMUXER.GAMMA_HIGH_GAIN =  100;
            C.SIMPLE_DEMUXER.GAMMA_LOW_GAIN  =    5;   % NOTE/POSSIBLE BUG: Uncertain which value is high-gain, and low-gain.
            
            obj.C = C;
        
            
            
            [obj.INPUTS_INFO_LIST,  obj.INPUTS_PDIDS_LIST]  = bicas.constants.produce_inputs_constants();
            [obj.OUTPUTS_INFO_LIST, obj.OUTPUTS_PDIDS_LIST] = bicas.constants.produce_outputs_constants(D);          
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
            error('BICAS:constants:Assertion', '"%s" is not a valid S/D mode ID', swModeId)
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
        function [outputsInfoList, eoutPdidList] = produce_outputs_constants(D)
            % TODO: Set SWD_LEVEL automatically?!
            
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
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s CWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s CWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'LFR L2s SWF science electric data in survey mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            % -------- TDS --------

            outputsInfoList{end+1} = [];            
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'TDS L2s CWF science electric data in low frequency mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            outputsInfoList{end+1} = [];
            outputsInfoList{end}.SWD_OUTPUT_FILE_IDENTIFIER = CLI_PARAMETER_SCI_NAME;
            outputsInfoList{end}.DATASET_ID                 = 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E';
            outputsInfoList{end}.SKELETON_VERSION_STR       = '02';
            outputsInfoList{end}.SWD_NAME                   =     'TDS L2s RSWF science electric data in low frequency mode';
            outputsInfoList{end}.SWD_DESCRIPTION            = 'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged';
            outputsInfoList{end}.SWD_LEVEL                  = 'L2S';
            outputsInfoList{end}.SWD_RELEASE_DATE           = D.INITIAL_RELEASE_DATE;
            outputsInfoList{end}.SWD_RELEASE_MODIFICATION   = D.INITIAL_RELEASE_MODIFICATION_STR;


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
