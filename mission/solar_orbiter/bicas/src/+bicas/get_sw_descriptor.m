% Return MATLAB structure that "exactly" corresponds to the S/W descriptor specified by the RCS ICD.
% Return result can be used for generating a string that can be printed.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created ~2016-06-01
%
%
% RETURN VALUE
% ============
% swDescriptor : See JSON_object_str for the exact format.
%
%
% IMPLEMENTATION NOTE
% ===================
% The values are derived entirely from BICAS constants. The structure is NOT directly incorporated in the BICAS
% constants/settings for extra flexibility. Reasons for NOT putting the S/W descriptor structure inside the settings
% class:
% (1) Some of the S/W descriptor variables have vague or misleading names ("name", "dataset versions", "dataset IDs")
%     which would (reasoably) have to be represented by MATLAB variables with the same names.
% (2) Some of the S/W descriptor variables are grouped in a way which
%     does not fit the rest of the code (e.g. modes[i].outputs.(output_XX).release in the S/W descriptor structure).
% (3) Some of the S/W descriptor values are really structure field NAMES, but would be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (4) The constants structure would become dependent on the format of the S/W descriptor structure.
%     The latter might change in the future, or be misunderstood in the present (i.e. be changed).
% (5) Some S/W descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (6) It is easier to add automatic checks on the S/W descriptor in the code that derives it.
%
function SwDescriptor = get_sw_descriptor(DataManager)
%
% PROPOSAL: Have this function add the prefix "input_" to all modes[].inputs[].input (in SWD)
% and only store the "CLI_PARAMETER_suffix" in the BICAS constants structure instead.
%   CON: The main function also uses the full CLI parameter. It too would have to add the prefix.
%        ==> The same "input_" prefix is specified in two places ==> Bad practice.
%
% PROPOSAL: Validation of the information.
%   PROPOSAL: Separate S/W descriptor validation code which always runs when BICAS is run.
%       PRO: Implementation (where the values are taken from) may change.
%   PROPOSAL: Use the JSON schema in the docs to verify the S/W descriptor automatically.
%       NOTE: Requires there to be code for validating.
%   PROPOSAL: All versions, CLI parameters, dates, dataset IDs on the right format.
%   PROPOSAL: No mode name, dataset ID doubles.
%   PROPOSAL: Check that there are no additional structure fields.
%   PROPOSAL: Check that paths, filenames are valid.
%   QUESTION: How handle overlap validation of BICAS constants, since all values come from there?
%   NOTE: Already implicitly checks that all the needed fields exist (since they are read here).
%   NOTE: Checks on the main constants structure will (can) only happen if this file is executed, not if 
%         the S/W as a whole is (by default).
%
global CONSTANTS SETTINGS

% Variable naming convention:
% ---------------------------
% Swd = S/W descriptor = The MATLAB structure which is used for producing the S/W descriptor (SWD) JSON object string.
% Its fields (field names) should NOT follow variable naming conventions since they determine the JSON object string
% which must follow the RCS ICD.

Swd = [];
Swd.identification.project     = SETTINGS.get('SWD_IDENTIFICATION.project');
Swd.identification.name        = SETTINGS.get('SWD_IDENTIFICATION.name');
Swd.identification.identifier  = SETTINGS.get('SWD_IDENTIFICATION.identifier');
Swd.identification.description = SETTINGS.get('SWD_IDENTIFICATION.description');
            
Swd.release.version            = SETTINGS.get('SWD_RELEASE.version');
Swd.release.date               = SETTINGS.get('SWD_RELEASE.date');
Swd.release.author             = SETTINGS.get('SWD_RELEASE.author');
Swd.release.contact            = SETTINGS.get('SWD_RELEASE.contact');
Swd.release.institute          = SETTINGS.get('SWD_RELEASE.institute');
Swd.release.modification       = SETTINGS.get('SWD_RELEASE.modification');

Swd.environment                = SETTINGS.get('SWD_ENVIRONMENT.executable');
Swd.modes = {};

for i = 1:length(CONSTANTS.SW_MODES_INFO_LIST)
    cliParameter = CONSTANTS.SW_MODES_INFO_LIST{i}.CLI_PARAMETER;
    
    ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(cliParameter);
    
    Swd.modes{end+1} = generate_sw_descriptor_mode(ExtendedSwModeInfo);
end


%===========================================================================================
% Validate S/W release version
% ----------------------------
% RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
% NOTE: It is hard to thoroughly follow the description, but the end result should be under
% release-->version-->pattern (not to be confused with release_dataset-->version--pattern).
%===========================================================================================
if isempty(regexp(Swd.release.version, '^(\d+\.)?(\d+\.)?(\d+)$', 'once'))
    error('BICAS:get_sw_descriptor:IllegalConfiguration', 'Illegal S/W descriptor release version "%s". This indicates a hard-coded configuration bug.', Swd.release.version)
end

SwDescriptor = Swd;    % Assign return value.

end



% Create data structure for a S/W mode corresponding to the information in the JSON S/W descriptor.
%
function SwdMode = generate_sw_descriptor_mode(ExtendedSwModeInfo)

global SETTINGS

SwdMode = struct;
SwdMode.name    = ExtendedSwModeInfo.CLI_PARAMETER;
SwdMode.purpose = ExtendedSwModeInfo.SWD_PURPOSE;

for OutputInfo = ExtendedSwModeInfo.inputs
    SwModeInputInfo = OutputInfo{1};
    
    SwdInput = struct;
    SwdInput.version    = SwModeInputInfo.SKELETON_VERSION_STR;
    SwdInput.identifier = SwModeInputInfo.DATASET_ID;
    
    SwdMode.inputs.(SwModeInputInfo.OPTION_HEADER_SH) = SwdInput;
end


for iOutput = 1:length(ExtendedSwModeInfo.outputs)
    OutputInfo    = ExtendedSwModeInfo.outputs{iOutput};

    [~, masterFilename] = bicas.get_master_CDF_path(OutputInfo.DATASET_ID, OutputInfo.SKELETON_VERSION_STR);

    SwdOutputInfo = struct;
    SwdOutputInfo.identifier           = OutputInfo.DATASET_ID;
    SwdOutputInfo.name                 = OutputInfo.SWD_NAME;
    SwdOutputInfo.description          = OutputInfo.SWD_DESCRIPTION;
    SwdOutputInfo.level                = OutputInfo.SWD_LEVEL;
    SwdOutputInfo.release.date         = OutputInfo.SWD_RELEASE_DATE;
    SwdOutputInfo.release.version      = OutputInfo.SKELETON_VERSION_STR;
    SwdOutputInfo.release.author       = SETTINGS.get('AUTHOR_NAME');
    SwdOutputInfo.release.contact      = SETTINGS.get('AUTHOR_EMAIL');
    SwdOutputInfo.release.institute    = SETTINGS.get('INSTITUTE');
    SwdOutputInfo.release.modification = OutputInfo.SWD_RELEASE_MODIFICATION;
    SwdOutputInfo.release.file         = masterFilename;
    
    % Validate output datasets release version
    % ----------------------------------------
    % RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
    % NOTE: It is hard to thoroughly follow the description, but the end result should be under
    % release_dataset-->version-->pattern (not to be confused with release-->version--pattern).
    if isempty(regexp(SwdOutputInfo.release.version, '^[0-9]{2}$', 'once'))
        error('BICAS:get_sw_descriptor:Assertion:IllegalConfiguration', ...
            'Illegal S/W descriptor output release version "%s". This indicates a hard-coded configuration bug.', SwdOutputInfo.release.version)
    end
        
    SwdMode.outputs.(OutputInfo.SWD_OUTPUT_FILE_IDENTIFIER) = SwdOutputInfo;
end



end
