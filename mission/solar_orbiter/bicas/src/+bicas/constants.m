% Constants - Singleton class for global constants used by BICAS.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software. Set up as a ~singleton handle class.
% Also contains validation code and functions for more convenient access.
%
%
% IMPORTANT NOTE
% ==============
% Some constants (1) correspond exactly to fields in the S/W (JSON) descriptor, and
% (2) are unlikely to be used for anything else. These are labeled with a prefix "SWD_". Other
% variables which do not have the prefix may also be used for the S/W descriptor too but they are
% probably more unambiguous in their meaning.
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
% PROPOSAL: Method for finding the path to the master CDF for a given output (instead of separate function).
%
% PROPOSAL: Include SW root path?! How?
%    PRO: Needs to be universally accessible.
%    CON: Not hardcoded. ==> Mixes code with manually set constants.
%    CON: No good MATLAB implementation of static variables.
%    QUESTION: Is there anything analogous? output dir?
%    PROPOSAL: Some functionality for setting "properties", in reality named ~global variables as key-value pairs. cf OVT.
%
% PROPOSAL: Merge with error_safe_constants?
%    CON: Would require to having only "error safe code".
%       PROPOSAL: Limit to code without branching (except lazy eval.) and without reading/calling
%                 external files so that one single clear+run will tell whether it works or not.
%          Ex: for loops are then safe.
%          CON: Can not use irf.log messages (for temp. data).
%       PROPOSAL: Have separate static method which is "error safe".
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
%
% PROPOSAL: Rename ".C" to ".general", ".misc"
% PROPOSAL: .inputs, .outputs to something more specific.
%   PROPOSAL: Something analogous to data_manager's "elementary input/output".
%   PROPOSAL: Something with CDF, datasets, BICAS_input/output, ...
%
% PROPOSAL: Change variable naming convention for C_* variables in general (not the inside constants.m).
%   Ex: C_sw_mode, C_input, C_output.
%   PROPOSAL: infoSwMode, infoInput, infoOutput
%   PROPOSAL: metadataSwMode, metadataInput, metadataOutput   (metaDataSwMode?)
%
% PROPOSAL: Convert constant cell arrays of structs to arrays of structs: S/W modes, CDF/PDID inputs/outputs.
%   PRO: Simplifies code that constructs cell arrays of the same struct field in multiple cell structs.
%
% PROPOSAL: assert_sw_mode_ID/assert_EI_PDID/assert_EO_PDID use lists of valid values.
%
% TODO: EI/EO --> EIn/EOut (must be done globally in other files?)
% TODO: Change more (all?) constants to upper case.
%###################################################################################################################

    properties(Access=public)
        C                  % For miscellaneous minor constants which still might require code to be initialized.
        BICAS_root_path
        %sw_descriptor
        inputs     % Information associated with input  datasets.
        outputs    % Information associated with output datasets.
        sw_modes   % Information associated with S/W modes datasets.
        EI_PDIDs
        EO_PDIDs
    end

    properties(Access=private)
        all_dataset_IDs    % Collect alla known dataset IDs. Useful for assertions.
    end
    
    %###################################################################################################################
    
    methods(Access=public)
        
        %=============
        % Constructor
        %=============
        function obj = constants(BICAS_root_path)            
            
            %--------------------------------------------------------------------------------
            % Common values. Only used INDIRECTLY to set the values of the "real" constants.
            %--------------------------------------------------------------------------------
            D = [];
            D.INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
            D.INITIAL_RELEASE_DATE = '2016-10-17';
            %D.SWD_OUTPUT_RELEASE_VERSION = '01';  % For the S/W descriptor output CDFs' release version. Unknown what a sensible value is.
            
            
            
            C = [];
            
            C.author_name         = 'Erik P G Johansson';
            C.author_email        = 'erik.johansson@irfu.se';
            C.institute           = 'IRF-U';
            C.master_CDFs_dir_rel = 'data';    % Location of master CDF files. Relative to the software directory structure root.
            
            % Value that shows up in EOut dataset GlobalAttributes.Calibration_version.
            % String value.
            C.Calibration_version = '0.1; Only proportionality constants; No voltage offset tables, no transfer functions';
        
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
            
            C.SWD_environment.executable = 'roc/bicas';
            
            % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
            C.stdout_prefix = 'STDOUT: ';
        
            % Parameters influencing how JSON objects are printed with function JSON_object_str.
            C.JSON_object_str = struct(...
                'indent_size',     4, ...
                'value_position', 15);
        
            % The epoch for ACQUISITION_TIME.
            % The time in UTC at which ACQUISITION_TIME is [0,0].
            % PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
            C.ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];

            % Define constants relating to LFR.
            % F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz.
            C.LFR = [];
            C.LFR.F0 = 24576;  % = 6 * 4096
            C.LFR.F1 =  4096;
            C.LFR.F2 =   256;
            C.LFR.F3 =    16;
        
            C.approximate_demuxer = struct(...
                'alpha',    1/17, ...
                'beta',        1, ...
                'gamma_hg',  100, ...
                'gamma_lg',    5);      % NOTE/POSSIBLE BUG: Uncertain which value is high-gain, and low-gain.
            
            C.INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID       = 0;   % Require input CDF Global Attribute "DATASET_ID"       to match the expected value.
            C.INPUT_CDF_ASSERTIONS.STRICT_Skeleton_version = 1;   % Require input CDF Global Attribute "Skeleton_version" to match the expected value.
            C.OUTPUT_CDF.SET_Test_id = 1;            % Set CDF GlobalAttribute "Test_id". ROC DFMD says that it should really be set by ROC.
            C.OUTPUT_CDF.Data_version = '01';        % Set CDF GlobalAttribute "Data_version". ROC DFMD says it should be updated in a way which can not be automatized?!!! Set here for now.
            
            C.PROCESSING.USE_AQUISITION_TIME_FOR_HK_INTERPOLATION = 1;
            
            % zVariables which are still empty after copying data into the master CDF assigned a correctly sized array with fill values.
            % This should only be necessary for S/W modes with incomplete processing.
            C.OUTPUT_CDF.EMPTY_ZVARIABLES_SET_TO_FILL = 1;
            
            C.IRF_LOG_LEVEL = 'notice';   % Log level for "irf.log".
            
            obj.C = C;
        
            
            
            [obj.inputs,  obj.EI_PDIDs]  = bicas.constants.produce_inputs_constants();
            [obj.outputs, obj.EO_PDIDs]  = bicas.constants.produce_outputs_constants(D);          
            obj.sw_modes = bicas.constants.produce_sw_modes_constants();
            
            obj.BICAS_root_path = BICAS_root_path;
            
            
            
            obj.all_dataset_IDs = unique(cellfun(@(s) ({s.dataset_ID}), [obj.outputs, obj.inputs])');

            
            obj.validate
        end
    end   % methods

    %###################################################################################################################

    methods(Access=public)
        
        function assert_dataset_ID(obj, dataset_ID)
        % Assert that argument is a valid dataset ID.
        
            if ~ismember(dataset_ID, obj.all_dataset_IDs)
                error('BICAS:constants:Assertion', '"%s" is not a valid dataset ID.', dataset_ID)
            end
        end
        
        function assert_sw_mode_ID(obj, sw_mode_ID)
            
            for i=1:length(obj.sw_modes)
                if strcmp(obj.sw_modes{i}.ID, sw_mode_ID)
                    return
                end
            end
            error('BICAS:constants:Assertion', '"%s" is not a valid S/D mode ID', sw_mode_ID)
        end

        function assert_EI_PDID(obj, EI_PDID)
            
            for i=1:length(obj.inputs)
                if strcmp(obj.inputs{i}.PDID, EI_PDID)
                    return
                end
            end
            error('BICAS:constants:Assertion', '"%s" is not a valid EI PDID', EI_PDID)
        end
%         
%         function assert_EO_PDID(obj, EO_PDID)
%             
%             for i=1:length(obj.outputs)
%                 if strcmp(obj.outputs{i}.PDID, EO_PDID)
%                     return
%                 end
%             end
%             error('BICAS:constants:Assertion', '"%s" is not a valid EO PDID', EO_PDID)
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
            for i = 1:length(obj.inputs)
                CLI_parameter = obj.inputs{i}.CLI_parameter;
                
                % NOTE: Implicitly checks that CLI_parameter does NOT begin with "--".
                disallowed_chars = setdiff(CLI_parameter, INPUT_CLI_PARAMETER_NAME_PERMITTED_CHARACTERS);
                if ~isempty(disallowed_chars)
                    error('BICAS:constants:Assertion:IllegalConfiguration', ...
                        'Constants value contains illegal character(s). This indicates a pure configuration bug (hard-coded).');
                end
            end
            
            bicas.utils.assert_strings_unique(obj.EI_PDIDs)
            bicas.utils.assert_strings_unique(obj.EO_PDIDs)            
            
            sw_mode_CLI_parameters = cellfun(@(s) ({s.CLI_parameter}), obj.sw_modes);
            sw_mode_IDs            = cellfun(@(s) ({s.ID           }), obj.sw_modes);
            bicas.utils.assert_strings_unique(sw_mode_CLI_parameters);
            bicas.utils.assert_strings_unique(sw_mode_IDs);
            
            % ASSERTION: CONSTANTS.sw_modes{i}.CLI_parameter matches validation regexp.
            for i = 1:length(obj.sw_modes)
                CLI_parameter = obj.sw_modes{i}.CLI_parameter;
                
                if isempty(regexp(CLI_parameter, SW_MODE_CLI_PARAMETER_REGEX, 'once'))
                    error('BICAS:constants:Assertion:IllegalConfiguration', ...
                        'Illegal S/W mode CLI parameter definition. This indicates a pure (hard-coded) configuration bug.');
                end
            end
            
            % NOTE: Check that combinations of dataset_ID and skeleton_version_str are unique.
            % Implemented by merging strings and checking for unique strings.
            % Is strictly speaking very slightly unsafe; could get false negatives.
            dataset_ID_version_list = cellfun( ...
                @(x) ({[x.dataset_ID, '_V', x.skeleton_version_str]}), ...
                [obj.outputs, obj.inputs]   );
            bicas.utils.assert_strings_unique(dataset_ID_version_list)
            
        end

    end   % methods(Access=private)
    
    %###################################################################################################################
    
    methods(Static, Access=private)
        
        
        function C_sw_modes = produce_sw_modes_constants()
            %===================================================================================================
            % Define the S/W modes which are visible to the "outside" of the software,
            % the modes which "officially" exist at any given time.
            %
            % Influences (at least) the required CLI arguments and the S/W descriptor.
            %
            % C_sw_modes : struct
            %    .CLI_parameter : Is used as CLI parameter to identify the S/W mode.
            %    .ID            : S/W mode ID. Used to identify the mode internally (in particular for hardcoded constants
            %                     in data_manager).
            %                     Has about the same purpose as CLI_parameter but is separate so that CLI_parameter
            %                     values/constants can be easily modified, whereas ID values are tied to hardcoded
            %                     constants in data_manager which are harder to modify.
            %
            %    .output_PDIDs : A cell array of PDIDs. Effectively an array of pointers to (1) the output constants, and (2)
            %                   indirectly to the input constants through data_manager.get_elementary_input_PDIDs.
            %===================================================================================================
            
            C_sw_modes = {};
            
            %=====
            % LFR 
            %=====
%             sw_mode = [];
%             sw_mode.CLI_parameter = 'LFR-SBM1-CWF-E';
%             sw_mode.ID            = 'LFR-SBM1-CWF-E_V01-V02';
%             sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR';            
%             sw_mode.output_PDIDs = {'L2S_LFR-SBM1-CWF-E_V02'};
%             C_sw_modes{end+1} = sw_mode;
%             
%             sw_mode = [];
%             sw_mode.CLI_parameter = 'LFR-SBM2-CWF-E';
%             sw_mode.ID            = 'LFR-SBM2-CWF-E_V01-V02';
%             sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR';            
%             sw_mode.output_PDIDs = {'L2S_LFR-SBM2-CWF-E_V02'};
%             C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SURV-CWF-E_V01-V02';
            sw_mode.ID            = 'LFR-SURV-CWF-E_V01-V02';
            sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR';
            sw_mode.input_PDIDs  = {'L2R_LFR-SURV-CWF_V01', 'HK_BIA_V01'};
            sw_mode.output_PDIDs = {'L2S_LFR-SURV-CWF-E_V02'};
            C_sw_modes{end+1} = sw_mode;
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'LFR-SURV-CWF-E_V02-V02';
            %sw_mode.ID            = 'LFR-SURV-CWF-E_V02-V02';
            %sw_mode.SWD_purpose = 'Generate CWF electric field data (potential difference) from LFR';            
            %sw_mode.input_PDIDs  = {'L2R_LFR-SURV-CWF_V02', 'HK_BIA_V01'};
            %sw_mode.output_PDIDs = {'L2S_LFR-SURV-CWF-E_V02'};
            %C_sw_modes{end+1} = sw_mode;
            
            sw_mode = [];
            sw_mode.CLI_parameter = 'LFR-SURV-SWF-E_V01-V02';
            sw_mode.ID            = 'LFR-SURV-SWF-E_V01-V02';
            sw_mode.SWD_purpose = 'Generate SWF electric (potential difference) data from LFR';            
            sw_mode.input_PDIDs  = {'L2R_LFR-SURV-SWF_V01', 'HK_BIA_V01'};
            sw_mode.output_PDIDs = {'L2S_LFR-SURV-SWF-E_V02'};
            C_sw_modes{end+1} = sw_mode;
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'LFR-SURV-SWF-E_V02-V02';
            %sw_mode.ID            = 'LFR-SURV-SWF-E_V02-V02';
            %sw_mode.SWD_purpose = 'Generate SWF electric (potential difference) data from LFR';
            %sw_mode.input_PDIDs  = {'L2R_LFR-SURV-SWF_V02', 'HK_BIA_V01'};
            %sw_mode.output_PDIDs = {'L2S_LFR-SURV-SWF-E_V02'};
            %C_sw_modes{end+1} = sw_mode;
            
            %=====
            % TDS
            %=====
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'TDS-LFM-CWF-E';
            %sw_mode.ID            = 'TDS-LFM-CWF-E_V01-V02';
            %sw_mode.SWD_purpose = 'Generate CWF electric (potential difference) data from TDS-LFM-CWF';
            %sw_mode.output_PDIDs = {'L2S_TDS-LFM-CWF-E_V02'};
            %C_sw_modes{end+1} = sw_mode;
            
            % NOTE: Accepts older/obsoleted V01 input data.
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'TDS-LFM-RSWF-E_V01-V02';
            %sw_mode.ID            = 'TDS-LFM-RSWF-E_V01-V02';
            %sw_mode.SWD_purpose = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V01';
            %sw_mode.output_PDIDs = {'L2S_TDS-LFM-RSWF-E_V02'};
            %C_sw_modes{end+1} = sw_mode;
            
            %sw_mode = [];
            %sw_mode.CLI_parameter = 'TDS-LFM-RSWF-E_V02-V02';
            %sw_mode.ID            = 'TDS-LFM-RSWF-E_V02-V02';
            %sw_mode.SWD_purpose = 'Generate RSWF electric (potential difference) data from TDS-LFM-RSWF V02';
            %sw_mode.output_PDIDs = {'L2S_TDS-LFM-RSWF-E_V02'};
            %C_sw_modes{end+1} = sw_mode;
            
        end
            
        
        
        %==========================================================
        % Produce constants for all possible INPUT datasets.
        % (independent of how they are associated with S/W modes).
        %==========================================================
        function [C_inputs, EI_PDIDs] = produce_inputs_constants
            % PROPOSAL: Put derivation of .PDID in nested init function.
            
            % NOTE: NESTED function
            function C_input = init_input_C_sci(dataset_ID, skeleton_version_str)
                C_input.CLI_parameter        = 'input_sci';
                C_input.dataset_ID           = dataset_ID;
                C_input.skeleton_version_str = skeleton_version_str;
            end

            C_inputs = {};
            
            % NOTE: CLI_parameter = CLI parameter MINUS flag prefix ("--"), i.e. e.g. "sci" instead of "--sci".
            
            %=========
            % BIAS HK
            %=========
            C_inputs{end+1} = [];
            C_inputs{end}.CLI_parameter        = 'input_hk';
            C_inputs{end}.dataset_ID           = 'ROC-SGSE_HK_RPW-BIA';
            C_inputs{end}.skeleton_version_str = '01';
            
            %=========
            % LFR SCI
            %=========            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '01');   % 1 snapshot/record
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM1-CWF', '02');   % 1   sample/record
            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '01');   % 1 snapshot/record
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SBM2-CWF', '02');   % 1   sample/record
            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '01');   % 1 snapshot/record
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-CWF', '02');   % 1   sample/record
            
            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '01');   % 1 snapshot/record
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-LFR-SURV-SWF', '02');   % 1 snapshot/record(!). Adds zVar SAMP_DTIME
            
            %=========
            % TDS SCI
            %=========
            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-CWF',  '01');  % 1   sample/record
            
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '01');  % 1   sample/record.  Adds two zVariables: SAMP_DTIME, SAMPS_PER_CH
            C_inputs{end+1} = init_input_C_sci('ROC-SGSE_L2R_RPW-TDS-LFM-RSWF', '02');  % 1 snapshot/record

            

            % Add one field ".PDID" to every struct above!
            % 
            % Put together PDIDs (used in data_manager).
            % See data_manager for definition.
            EI_PDIDs = {};
            for i = 1:length(C_inputs)
                C_inputs{i}.PDID = bicas.constants.construct_PDID(C_inputs{i}.dataset_ID, C_inputs{i}.skeleton_version_str);
                EI_PDIDs{i} = C_inputs{i}.PDID;
            end
        end



        %==========================================================
        % Produce constants for all possible OUTPUT datasets
        % (independent of how they are associated with S/W modes).
        %==========================================================
        function [C_outputs, EO_PDIDs] = produce_outputs_constants(D)
            % TODO: Set SWD_level automatically?!
            
            CLI_PARAMETER_SCI_NAME = 'output_sci';
            
            C_outputs = {};
            
            % -------- LFR --------
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 1, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in selective burst mode 2, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'LFR L2s CWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'LFR L2s SWF science electric data in survey mode';
            C_outputs{end}.SWD_description             = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            % -------- TDS --------

            C_outputs{end+1} = [];            
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'TDS L2s CWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s CWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;
            
            C_outputs{end+1} = [];
            C_outputs{end}.JSON_output_file_identifier = CLI_PARAMETER_SCI_NAME;
            C_outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E';
            C_outputs{end}.skeleton_version_str        = '02';
            C_outputs{end}.SWD_name                    =     'TDS L2s RSWF science electric data in low frequency mode';
            C_outputs{end}.SWD_description             = 'RPW TDS L2s RSWF science electric (potential difference) data in low frequency mode, time-tagged';
            C_outputs{end}.SWD_level                   = 'L2S';
            C_outputs{end}.SWD_release_date            = D.INITIAL_RELEASE_DATE;
            C_outputs{end}.SWD_release_modification    = D.INITIAL_RELEASE_MODIFICATION_STR;
            %C_outputs{end}.SWD_release_version         = D.SWD_OUTPUT_RELEASE_VERSION;


            % Add one field ".PDID" to every struct above!
            % 
            % Put together PDIDs (used in data_manager).
            % See data_manager for definition.
            EO_PDIDs = {};
            for i = 1:length(C_outputs)
                C_outputs{i}.PDID = bicas.constants.construct_PDID(C_outputs{i}.dataset_ID, C_outputs{i}.skeleton_version_str);
                EO_PDIDs{i} = C_outputs{i}.PDID;
            end
        end
        
        
        
        function PDID = construct_PDID(dataset_ID, skeleton_version_str)
        % Construct a (shorter) PDID from a dataset ID and skeleton version).
        
            %PDID = [dataset_ID, '_V', skeleton_version_str];
            
            dataset_ID_shortened = regexprep(dataset_ID,           '^ROC-SGSE_', '',  'once');
            dataset_ID_shortened = regexprep(dataset_ID_shortened, '_RPW-',      '_', 'once');                
            PDID = [dataset_ID_shortened, '_V', skeleton_version_str];
        end

    end % methods(Static, Access=private)
    
end   % classdef
