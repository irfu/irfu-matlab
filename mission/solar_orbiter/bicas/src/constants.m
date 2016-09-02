% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software. Set up as a ~singleton handle class.
% Also contains validation code and functions for more convenient access.
%
% IMPORTANT NOTE: Some constants (1) correspond exactly to fields in the S/W (JSON) descriptor, and
% (2) are unlikely to be used for anything else. These are labeled with a prefix "SWD_". Other
% variables which do not have the prefix may also be used for the S/W descriptor too but they are
% probably more unambiguous in their meaning.
%
% IMPLEMENTATION NOTE: Reasons for using a singleton class (instead of static methods:
% 1) Can use properties/instance variables for "caching" values. Persistent variables are bad for testing.
% NOTE: There are no proper static variables in MATLAB.
% 2) Can split up (structure, organize) configuration and validation code in methods.
% 3) The constructor can be used as initialization code which must be run before using the class/constants.
%
% IMPLEMENTATION NOTE: There is other code which builds a structure corresponding to the S/W descriptor
% from the constants structure here.
% Reasons for NOT putting the S/W descriptor structure inside the constants structure:
% (1) Some of the S/W descriptor variables have vague or misleading names ("name", "dataset versions", "dataset IDs")
%     which would be represented by MATLAB variables by the same names.
% (2) Some of the S/W descriptor variables are grouped in a way which
%     does not fit the rest of the code (modes[].outputs.output_XX.release in the S/W descriptor structure).
% (2) Some of the S/W descriptor values are really structure field NAMES, but would be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (3) The constants structure would become dependent on the format of the S/W descriptor structure.
%     The latter might change in the future, or be misunderstood in the precent (i.e. be changed).
% (4) Some S/W descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (5) It is easier to add automatic checks on the S/W descriptor in the code that derives it.
%
% NOTE: sw_root_dir is not strictly a constant.
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
% PROPOSAL: Method for finding the path to the master CDF for a given output.
%
% PROPOSAL: Include SW root path?! How?
%    PRO: Needs to be universally accessible.
%    CON: Not hardcoded. ==> Mixes code with manually set constants.
%    CON: No good MATLAB implementation of static variables.
%    QUESTION: Is there anything analogous? output dir?
%
% PROPOSAL: Merge with init global constants?
%    CON: Calling would slow down code. Many functions use ERROR_CODES.
%       PROPOSAL: Use one big global variable instead (initialized here).
%    CON: Would require to having only "error safe code".
%       PROPOSAL: Limit to code without branching (except lazy eval.) and without reading/calling
%                 external files so that one single clear+run will tell whether it works or not.
%          Ex: for loops are then safe.
%          CON: Can not use irf.log messages (for temp. data).
%       PROPOSAL: Have separate static method which is "error safe".
% PROPOSAL: Add S/W descriptor to this structure by calling get_sw_descriptor?
%    PRO: Can force validation to take place in get_sw_descriptor.
%
% PROPOSAL: More validation
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



    % NOTE: Should be private when not debugging.
    properties(Access=private)
        inputs
        outputs
        root_dir_path
    end
    
    %###################################################################################################################
    
    properties(Access=public)
        C                  % For miscellaneous minor constants which still might require code to be initialized.
        sw_modes
    end

    %###################################################################################################################
    
    methods(Access=public)
        
        %=============
        % Constructor
        %=============
        function obj = constants()            
            
            %--------------------------------------------------------------------------------
            % Common values. Only used INDIRECTLY to set the values of the "real" constants.
            %--------------------------------------------------------------------------------
            D = [];
            D.INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            D.INITIAL_RELEASE_DATE = '2016-09-01';
            %D.SWD_OUTPUT_RELEASE_VERSION = '01';  % For the S/W descriptor output CDFs' release version. Unknown what a sensible value is.
            
            
            
            C = [];
            
            C.author_name         = 'Erik P G Johansson';
            C.author_email        = 'erik.johansson@irfu.se';
            C.institute           = 'IRF-U';
            C.master_CDFs_dir_rel = 'data';    % Location of master CDF files. Relative to the software directory structure root.
        
            %irf.log('w', 'Using temporary S/W name in constants.')
            C.SWD_identification.project     = 'ROC-SGSE';
            C.SWD_identification.name        = 'BICAS';
            C.SWD_identification.identifier  = 'ROC-SGSE-BICAS';
            C.SWD_identification.description = 'BIAS Calibration Software (BICAS) which derives the BIAS L2S input signals (plus some) from the BIAS L2R output signals.';
            
            % Refers to the S/W descriptor release data for the entire software (not specific outputs).
            C.SWD_release.version      = '0.1.0';
            C.SWD_release.date         = D.INITIAL_RELEASE_DATE;
            C.SWD_release.author       = C.author_name;
            C.SWD_release.contact      = C.author_email;
            C.SWD_release.institute    = C.institute;
            C.SWD_release.modification = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            C.SWD_environment.executable = 'roc/bicas';   % Temporary SW (file) name
            
            % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
            C.stdout_prefix = 'STDOUT: ';
        
            % Parameters influencing how JSON objects are printed with function JSON_object_str.
            C.JSON_object_str = struct(...
                'indent_size',     4, ...
                'value_position', 15);
        
            % Define constants relating to LFR.
            % F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz.
            C.LFR = [];
            C.LFR.F0 = 24576;  % = 6 * 4096
            C.LFR.F1 =  4096;
            C.LFR.F2 =   256;
            C.LFR.F3 =    16;
        
            C.approximate_demuxer = struct(...
                'alpha',    1/17, ...
                'beta',       1, ...
                'gamma_hg',   5, ...
                'gamma_lg', 100);
            
            obj.C = C;
        
            
            
            obj.inputs   = constants.produce_inputs_constants();
            obj.outputs  = constants.produce_outputs_constants(D);          
            obj.sw_modes = constants.produce_sw_modes_constants();
            
            
            
            obj.validate
        end
    end   % methods
    
    
    %###################################################################################################################
    
    methods(Access=public)
        
        %===================================================================================================
        
        %=========================================================================
        % Effectively a variable that can be written to once, and then only read.
        %=========================================================================
        function varargout = SW_root_dir(obj, varargin)
            % PROPOSAL: Change to plain public variable/property?
            global ERROR_CODES
            
            if nargout == 0 && length(varargin) == 1 && isempty(obj.root_dir_path)
                obj.root_dir_path = varargin{1};
            elseif length(varargin) == 0 && ~isempty(obj.root_dir_path)  % NOTE: Does not check for nargout intentionally.
                varargout{1} = obj.root_dir_path;
            else
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Trying to set already set constant, or reading unset constant. Pure code bug.')
            end
        end

        
        
        %================================================================================================================
        % Return constants structure for a specific S/W mode.
        %        
        % This structure is automatically put together from other structures (constants) to avoid having to define too
        % many redundant constants.
        % NOTE: The function takes the CLI_parameter as parameter, not the ID.
        %================================================================================================================
        function C_sw_mode = get_C_sw_mode_full(obj, CLI_parameter)
            global ERROR_CODES
            
            C_sw_mode = select_structs(obj.sw_modes, 'CLI_parameter', {CLI_parameter});            
            C_sw_mode = C_sw_mode{1};
            
            %=================================================
            % Collect all elementary input process data types
            %=================================================
            input_process_data_types = {};            
            for i = 1:length(C_sw_mode.output_process_data_types)
                %==============================================================================================
                % Find all elementary input process data types for a given elementary output process data type
                %==============================================================================================
                temp = data_manager.get_elementary_input_process_data_types(...
                    C_sw_mode.output_process_data_types{i}, C_sw_mode.ID);
                
                input_process_data_types = [input_process_data_types, temp];
            end
            
            try
                C_sw_mode.inputs  = select_structs(obj.inputs,  'process_data_type', input_process_data_types);
            catch exception
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Can not identify all input process data types associated with mode/CLI parameter "%s".', CLI_parameter)
            end
            try
                C_sw_mode.outputs = select_structs(obj.outputs, 'process_data_type', C_sw_mode.output_process_data_types);
            catch exception
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Can not identify all output process data types associated with mode/CLI parameter "%s".', CLI_parameter)
            end
        end

    end   % methods
    
    %###################################################################################################################
    
    methods(Access=private)

        % Any code for double-checking the validity of hardcoded constants.
        function validate(obj)

            global ERROR_CODES            
            C = obj.C;



            %==========================
            % Iterate over input types
            %==========================
            % The RCS ICD, iss2rev2, section 3.2.3 only permits these characters (and only lowercase).
            INPUT_CLI_PARAMETER_NAME_PERMITTED_CHARACTERS = 'abcdefghijklmnopqrstuvxyz0123456789_';
            for i = 1:length(obj.inputs)
                CLI_parameter = obj.inputs{i}.CLI_parameter;
                
                % NOTE: Implicitly checks that CLI_parameter does NOT begin with "--".
                disallowed_chars = setdiff(CLI_parameter, INPUT_CLI_PARAMETER_NAME_PERMITTED_CHARACTERS);
                if ~isempty(disallowed_chars)
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Constants value contains illegal character(s). This indicates a pure configuration bug (hard-coded).');
                end
            end            
            
            %========================
            % Iterate over S/W modes
            %========================
            sw_mode_CLI_parameters = {};
            sw_mode_IDs            = {};
            % The RCS ICD, iss2rev2, section 5.3 seems (ambiguous) to imply this regex.
            SW_MODE_CLI_PARAMETER_REGEX = '^[A-Za-z][\w-]+$';   % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.
            for i = 1:length(obj.sw_modes)
                sw_mode_CLI_parameter = obj.sw_modes{i}.CLI_parameter;
                
                sw_mode_CLI_parameters{end+1} = sw_mode_CLI_parameter;
                sw_mode_IDs{end+1}            = obj.sw_modes{i}.ID;
                C_sw_mode = obj.get_C_sw_mode_full(sw_mode_CLI_parameter);
                
                if ~length(regexp(sw_mode_CLI_parameter, SW_MODE_CLI_PARAMETER_REGEX))
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal S/W mode CLI parameter definition. This indicates a pure configuration bug (hard-coded).');
                end
                
                % Iterate over inputs
                inputs_CLI_parameters = {};
                for j = 1:length(C_sw_mode.inputs)
                    inputs_CLI_parameters{end+1} = C_sw_mode.inputs{j}.CLI_parameter;
                end
                assert_strings_unique(inputs_CLI_parameters)
                
                % Iterate over outputs
                outputs_JSON_output_file_identifiers = {};
                for j = 1:length(C_sw_mode.outputs)
                    outputs_JSON_output_file_identifiers{end+1} = C_sw_mode.outputs{j}.JSON_output_file_identifier;
                end
                assert_strings_unique(outputs_JSON_output_file_identifiers)
            end
            assert_strings_unique(sw_mode_CLI_parameters)
            assert_strings_unique(sw_mode_IDs)

            
            
            % NOTE: Check that combinations of dataset_ID and skeleton_version_str are unique.
            % Implemented by merging strings and checking for unique strings.
            % Is strictly speaking very slightly unsafe; could get false negatives.
            dataset_ID_version_list = cellfun( ...
                @(x) ({[x.dataset_ID, '_V', x.skeleton_version_str]}), ...
                [obj.outputs, obj.inputs]   );
            assert_strings_unique(dataset_ID_version_list)           
        end

    end   % methods
    
    %###################################################################################################################
    
    methods(Static, Access=private)
        
        
        function C_sw_modes = produce_sw_modes_constants()
            %---------------------------------------------------------------------------------------------------
            % Define the S/W modes which are visible to the "outside" of the software,
            % the modes which "officially" exist at any given time.
            %
            % Influences (at least) the required CLI arguments and the S/W descriptor.
            %---------------------------------------------------------------------------------------------------
            % .CLI_parameter : Is used as CLI parameter to identify the S/W mode.
            % .ID            : S/W mode ID. Used to identify the mode internally (in particular for hardcoded constants
            %                  in data_manager).
            %                  Has about the same purpose as CLI_parameter but is separate so that CLI_parameter
            %                  values/constants can be easily modified, whereas ID values are tied to hardcoded
            %                  constants in data_manager which are harder to modify.
            %
            % .output_process_data_types : Effectively an array of pointers to (1) the output constants, and (2)
            % indirectly to the input constants through data_manager.get_elementary_input_process_data_types.
            
            C_sw_modes = {};
            
            % -------- LFR --------
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SURV-CWF-E';
            sw_mode.ID            = 'LFR-SURV-CWF-E_V01-V01';
            sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR';            
            sw_mode.output_process_data_types = {'L2S_LFR-SURV-CWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SURV-CWF-E___EXPERIMENTAL';
            sw_mode.ID            = 'LFR-SURV-CWF-E_V01-V01___EXPERIMENTAL';
            sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR - EXPERIMENTAL';            
            sw_mode.output_process_data_types = {'L2S_LFR-SURV-CWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SURV-SWF-E';
            sw_mode.ID            = 'LFR-SURV-SWF-E_V01-V01';
            sw_mode.SWD_purpose = 'Generate SWF electric (potential difference) data from LFR';            
            sw_mode.output_process_data_types = {'L2S_LFR-SURV-SWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            % -------- TDS --------
            sw_mode = [];
            sw_mode.CLI_parameter = 'TDS-LFM-CWF-E';
            sw_mode.ID            = 'TDS-LFM-CWF-E_V01-V01';
            sw_mode.SWD_purpose = 'Generate CWF electric (potential difference) data from TDS-LFM-CWF';
            sw_mode.output_process_data_types = {'L2S_TDS-LFM-CWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'TDS-LFM-RSWF-E_V01-V01';
            sw_mode.ID            = 'TDS-LFM-RSWF-E_V01-V01';
            sw_mode.SWD_purpose = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V01';
            sw_mode.output_process_data_types = {'L2S_TDS-LFM-RSWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'TDS-LFM-RSWF-E_V02-V01';
            sw_mode.ID            = 'TDS-LFM-RSWF-E_V02-V01';
            sw_mode.SWD_purpose = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V02';
            sw_mode.output_process_data_types = {'L2S_TDS-LFM-RSWF-E_V01'};
            C_sw_modes{end+1} = sw_mode;
            
            % -------- TEST --------
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'TEST-MODE';
            %sw_mode.SWD_purpose = 'Test mode. This mode is never supposed to be seen outside of development.';            
            %sw_mode.output_process_data_types = {'L2S_LFR-SURV-SWF-E_V01', 'L2S_TEST_V99'};
            %C_sw_modes{end+1} = sw_mode;
        end
            
        
        %==========================================================
        % Produce constants for all possible INPUT datasets.
        % (independent of how they are associated with S/W modes).
        %==========================================================
        function C_inputs = produce_inputs_constants
            CLI_PARAMETER_SCI_NAME = 'input_sci';
            
            C_inputs = {};
            
            % NOTE: CLI_parameter = CLI parameter MINUS flag prefix ("--").           
            
            %=========
            % BIAS HK
            %=========
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = 'input_hk';
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_HK_RPW-BIA';
            C_inputs{end}.skeleton_version_str = '01';
            
            %=========
            % LFR SCI
            %=========
            % 1 snapshot/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-LFR-SBM1-CWF';
            C_inputs{end}.skeleton_version_str = '01';
            
            % 1 snapshot/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-LFR-SBM2-CWF';
            C_inputs{end}.skeleton_version_str = '01';
            
            % 1 snapshot/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF';
            C_inputs{end}.skeleton_version_str = '01';
            
            % 1 snapshot/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF';
            C_inputs{end}.skeleton_version_str = '01';

            %=========
            % TDS SCI
            %=========
            % 1 sample/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-TDS-LFM-CWF';
            C_inputs{end}.skeleton_version_str = '01';

            % 1 sample/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-TDS-LFM-RSWF';
            C_inputs{end}.skeleton_version_str = '01';
            
            % 1 snapshot/record
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_RPW-TDS-LFM-RSWF';
            C_inputs{end}.skeleton_version_str = '02';
            
            %======
            % TEST
            %======
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter        = 'input_sci2';
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_L2R_TEST';
            C_inputs{end}.skeleton_version_str = '99';
            
            for i = 1:length(C_inputs)
                C_inputs{i}.process_data_type = constants.construct_process_data_type(C_inputs{i}.dataset_ID, C_inputs{i}.skeleton_version_str);
                %C_inputs{i}.process_data_type = [C_inputs{i}.dataset_ID, '_V', C_inputs{i}.skeleton_version_str];
            end
        end
        
        
        
        %==========================================================
        % Produce constants for all possible OUTPUT datasets
        % (independent of how they are associated with S/W modes).
        %==========================================================
        function C_outputs = produce_outputs_constants(D)
            % TODO: Set SWD_level automatically?!
            
            CLI_PARAMETER_SCI_NAME = 'output_sci';
            
            C_outputs = {};
            
            % -------- LFR --------
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'LFR L2s SWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 1, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            % -------- TDS --------

            C_outputs{end+1} = [];            
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'TDS L2s CWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E';
            C_outputs{end}.skeleton_version_str        = '01';
            C_outputs{end}.SWD_name                    =     'TDS L2s RSWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            % -------- TEST --------
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = 'output_test2';
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_TEST';
            C_outputs{end}.skeleton_version_str        = '99';
            C_outputs{end}.SWD_name                    = 'Test form of output.';
            C_outputs{end}.SWD_description             = 'Test form of output. This never supposed to be seen outside of development.';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            % Put together process data types (used in data_manager).
            % See data_manager for definition.
            for i = 1:length(C_outputs)
                %C_outputs{i}.process_data_type = [C_outputs{i}.dataset_ID, '_V', C_outputs{i}.skeleton_version_str];
                C_outputs{i}.process_data_type = constants.construct_process_data_type(C_outputs{i}.dataset_ID, C_outputs{i}.skeleton_version_str);
            end
        end
        
        
        
        function process_data_type = construct_process_data_type(dataset_ID, skeleton_version_str)
            %process_data_type = [dataset_ID, '_V', skeleton_version_str];
            
            dataset_ID_shortened = regexprep(dataset_ID,           '^ROC-SGSE_', '',  'once');
            dataset_ID_shortened = regexprep(dataset_ID_shortened, '_RPW-',      '_', 'once');                
            process_data_type = [dataset_ID_shortened, '_V', skeleton_version_str];
            %disp(process_data_type)   % DEBUG
        end

    end % methods
    
end   % classdef
