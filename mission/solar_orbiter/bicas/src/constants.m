% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software. Set up as a ~singleton handle class.
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
% Reasons for not putting the S/W descriptor structure inside the constants structure:
% (1) Some of the S/W descriptor variables have vague or misleading names ("name", "dataset versions", "dataset IDs")
% (2) Some of the S/W descriptor variables are grouped in a way which
%     does not fit the rest of the code (modes[].outputs.output_XX.release in the S/W descriptor structure).
% (2) Some of the S/W descriptor values are really structure field NAMES, but should be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (3) The constants structure would become dependent on the format of the S/W descriptor structure.
%     The latter might change in the future.
% (4) Some S/W descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (5) It is easier to add automatic checks on the S/W descriptor in the code that derives it.
%
% NOTE: sw_root_dir is not strictly a constant.
%
classdef constants < handle
%
% PROPOSAL: Redefine file into something that makes "decisions"?
% PROPOSAL: Method for finding the path to the master cdf for a given output.
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
%   PROPOSAL: Check that master cdfs exist, that paths exist.
%       CON: Makes sense to make that kind of check here?!
%   PROPOSAL: Check that data types are unique.
%       NOTE: Requires access to the lists.
%
% PROPOSAL: Have string identifiers for input/output types identical to what "data_manager" uses?
%    PROPOSAL: Use containers.Map instead of cell arrays.
%       CON: Using dataset IDs as keys creates another double use of that value.
%       PRO: Makes it more convenient to internally have multiple input/output dataset versions.
%


    % NOTE: Should be private when not debugging.
    properties(Access=private)
    %properties     
        C
        inputs
        outputs
        root_dir_path
    end
    
    properties(Access=public,Constant)
        
        % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
        stdout_prefix = 'STDOUT: ';
        
        % Parameters influencing how JSON objects are printed with function JSON_object_str.
        JSON_object_str = struct(...
            'indent_size',          4, ...
            'fill_str_max_length', 13);
        
        approximate_demuxer = struct(...
            'alpha',    1/17, ...
            'beta',       1, ...
            'gamma_hg',   5, ...
            'gamma_lg', 100);
    end

    methods(Access=public)
        
        % Constructor
        function obj = constants()            
            % PROPOSAL: Other name for sw_mode.CLI_parameter which implies generic string identifier (which is
            %           used to derive a CLI parameter).
            
            %--------------------------------------------------------------------------------
            % Common values. Only used INDIRECTLY to set the values of the "real" constants.
            %--------------------------------------------------------------------------------
            D = [];
            D.INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            D.INITIAL_RELEASE_DATE = '2016-07-21';
                        
            
            
            obj.inputs = constants.produce_inputs_constants();
            obj.outputs = constants.produce_outputs_constants(D);          
            
            

            C = [];
            C.author_name = 'Erik P G Johansson';
            C.author_email = 'erik.johansson@irfu.se';
            C.institute = 'IRF-U';
            C.master_cdfs_dir_rel = 'data';    % Location of master CDF files. Relative to the software directory structure root.
            
            
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
            
            C.sw_modes = {};
            
            
            
            %---------------------------------------------------------------------------------------------------
            % Define the software modes which are presented "outside" the software.
            %
            % Influences (at least) the required CLI arguments and the S/W descriptor.
            %---------------------------------------------------------------------------------------------------
            % NOTE: sw_mode.CLI_parameter : Is used as CLI parameter, but also to identify the mode.
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-CWF-E';
            sw_mode.SWD_purpose = 'Generate continuous waveform electric field data (potential difference) from LFR';            
            sw_mode.inputs  = select_structs(obj.inputs,  'dataset_ID', {'ROC-SGSE_L2R_RPW-LFR-SURV-CWF', 'ROC-SGSE_HK_RPW-BIA'});
            sw_mode.outputs = select_structs(obj.outputs, 'dataset_ID', {'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E'});
            C.sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SWF-E';
            sw_mode.SWD_purpose = 'Generate snapshow waveform electric (potential difference) data from LFR';            
            sw_mode.inputs  = select_structs(obj.inputs,  'dataset_ID', {'ROC-SGSE_L2R_RPW-LFR-SURV-SWF', 'ROC-SGSE_HK_RPW-BIA'});
            sw_mode.outputs = select_structs(obj.outputs, 'dataset_ID', {'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E'});            
            C.sw_modes{end+1} = sw_mode;

            
            % Define constants relating to LFR.
            % F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz.
            C.LFR = [];
            C.LFR.F0 = 24576;   % = 6 * 4096
            C.LFR.F1 = 4096;
            C.LFR.F2 = 256;
            C.LFR.F3 = 16;
            
            
            
            obj.C = C;
            
            obj.validate
        end
        
        %===================================================================================================
        
        % Effectively a variable that can be written to once, and then only read.
        % Change to plain public variable/property?
        function varargout = SW_root_dir(obj, varargin)
            global ERROR_CODES
            
            if nargout == 0 && length(varargin) == 1 && isempty(obj.root_dir_path)
                obj.root_dir_path = varargin{1};
            elseif length(varargin) == 0 && ~isempty(obj.root_dir_path)  % NOTE: Does not check for nargout intentionally.
                varargout{1} = obj.root_dir_path;
            else
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Trying to set already set constant, or reading unset constant. Pure code bug.')
            end
        end
        
        %===================================================================================================
        
        function C = get_general(obj)
            C = obj.C;
        end
        
        % Return info on a subset of possible INPUT data sets.
        function inputs = get_cdf_inputs_constants(obj, dataset_IDs)           
            inputs = select_structs(obj.inputs, 'dataset_ID', dataset_IDs);            
        end
        
        % Return info on a subset of possible OUTPUT data sets.
        function outputs = get_cdf_outputs_constants(obj, dataset_IDs)           
            outputs = select_structs(obj.outputs, 'dataset_ID', dataset_IDs);            
        end
        
    end   % methods
    
    
    
    methods(Static, Access=private)
        
        % Produce constants for all possible INPUT data sets.
        function C_inputs = produce_inputs_constants
            CLI_PARAMETER_SCI_NAME = 'input_sci';
            
            C_inputs = {};
            
            % NOTE: CLI_parameter_name = CLI parameter MINUS flag prefix ("--").           
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = 'input_hk';
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_HK_RPW-BIA';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SBM1-CWF';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SBM2-CWF';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-TDS-LFM-RSWF';
            C_inputs{end}.dataset_version_str = '01';
            
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter_name  = CLI_PARAMETER_SCI_NAME;
            C_inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-TDS-LFM-CWF';
            C_inputs{end}.dataset_version_str = '01';
            
        end
        
        % Produce constants for all possible OUTPUT data sets.
        function C_outputs = produce_outputs_constants(D)
            % PROPOSAL: Abolish master_cdf_filename? Return the value from a function?
            %   CON: The value appears in the S/W descriptor and is hence "registered". There could
            %   be some value in having an actual table.
            %
            % TODO: Set master_cdf_filename automatically.
            % TODO: Set SWD_level automatically?!
            
            CLI_PARAMETER_SCI_NAME = 'output_sci';
            
            C_outputs = {};
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01.cdf';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01.cdf';
            C_outputs{end}.SWD_name                    =     'LFR L2s SWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = '';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 1, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = '';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;

            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = '';
            C_outputs{end}.SWD_name                    =     'TDS L2s RSWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            C_outputs{end+1} = [];
            
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E';
            C_outputs{end}.dataset_version_str         = '01';
            %C_outputs{end}.master_cdf_filename         = '';
            C_outputs{end}.SWD_name                    =     'TDS L2s CWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            
            for i=1:length(C_outputs)
                C_outputs{i}.master_cdf_filename = [C_outputs{i}.dataset_ID, '_V', C_outputs{end}.dataset_version_str, '.cdf'];
                %C_outputs{i}.SWD_level = 'L2S';
            end
        end
        
    end % methods
    
    
    
    methods(Access=private)

        % Any code for double-checking the validity of hardcoded constants.
        function validate(obj)

            global ERROR_CODES
            
            C = obj.C;

            for sw_mode = C.sw_modes
                % The RCS ICD, iss2rev2, section 3.2.3 only permits these characters (and only lowercase).
                PERMITTED_CHARACTERS = 'abcdefghijklmnopqrstuvxyz0123456789_';
                
                % NOTE: Implicitly checks that CLI_parameter_name does not begin with "--".
                for input = sw_mode{1}.inputs
                    disallowed_chars = setdiff(input{1}.CLI_parameter_name, PERMITTED_CHARACTERS);
                    if ~isempty(disallowed_chars)
                        errorp(ERROR_CODES.ASSERTION_ERROR, 'Constants value contains illegal characters. This indicates a pure bug.');
                    end
                end
            end
            
            % NOTE: Can not handle multiple version of same dataset ID. Should validate combinations
            % thereof.
            dataset_IDs = cellfun(@(x) ({x.dataset_ID}), {obj.outputs{:}, obj.inputs{:}});                        
            validate_strings_unique(dataset_IDs)           
        end

    end   % methods    
    
end   % classdef
