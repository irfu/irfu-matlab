% Return MATLAB structure that "exactly" corresponds to the S/W descriptor specified by the RCS ICD.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created ~2016-06-01
%
% The values are derived entirely from BICAS constants. The structure is NOT directly incorporated in the BICAS
% constants for extra flexibility.
%
function swDescriptor = get_sw_descriptor(DataManager)
%
% PROPOSAL: Have this function add the prefix "input_" to all modes[].inputs[].input (in SWD)
% and only store the "CLI_parameter_suffix" in the BICAS constants structure instead.
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
global CONSTANTS

% SWD = The structure which is used for producing the S/W descriptor (SWD) JSON object string. Its fields (field names)
% should NOT follow variable naming conventions since they influence the JSON object string.
swd = [];
swd.identification = CONSTANTS.C.SWD_IDENTIFICATION;
swd.release        = CONSTANTS.C.SWD_RELEASE;
swd.environment    = CONSTANTS.C.SWD_ENVIRONMENT;
swd.modes = {};

for i = 1:length(CONSTANTS.SW_MODES_INFO_LIST)
    cliParameter = CONSTANTS.SW_MODES_INFO_LIST{i}.CLI_parameter;
    
    ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(cliParameter);
    
    swd.modes{end+1} = generate_sw_descriptor_mode(CONSTANTS.C, ExtendedSwModeInfo);
end


%===========================================================================================
% Validate S/W release version
% ----------------------------
% RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
% NOTE: It is hard to thoroughly follow the description, but the end result should be under
% release-->version-->pattern (not to be confused with release_dataset-->version--pattern).
%===========================================================================================
if isempty(regexp(swd.release.version, '^(\d+\.)?(\d+\.)?(\d+)$', 'once'))
    error('BICAS:get_sw_descriptor:IllegalConfiguration', 'Illegal S/W descriptor release version "%s". This indicates a hard-coded configuration bug.', swd.release.version)
end

swDescriptor = swd;

end



%===================================================================================================
% Create data structure for a S/W mode corresponding to the information in the JSON S/W descriptor.
%===================================================================================================
function swdMode = generate_sw_descriptor_mode(C, ExtendedSwModeInfo)
%
% Variable naming convention:
%    SWD = S/W descriptor

swdMode = [];
swdMode.name    = ExtendedSwModeInfo.CLI_parameter;
swdMode.purpose = ExtendedSwModeInfo.SWD_purpose;

for x = ExtendedSwModeInfo.inputs
    swModeInput = x{1};
    
    swdInput = [];    
    swdInput.version    = swModeInput.skeleton_version_str;
    swdInput.identifier = swModeInput.dataset_ID;
    
    swdMode.inputs.(swModeInput.CLI_parameter) = swdInput;
end

for x = ExtendedSwModeInfo.outputs
    swModeOutput = x{1};
    swdOutput = [];

    [~, masterFilename] = bicas.get_master_CDF_path(swModeOutput.dataset_ID, swModeOutput.skeleton_version_str);

    swdOutput.identifier  = swModeOutput.dataset_ID;
    swdOutput.name        = swModeOutput.SWD_name;
    swdOutput.description = swModeOutput.SWD_description;
    swdOutput.level       = swModeOutput.SWD_level;
    swdOutput.release.date         = swModeOutput.SWD_release_date;
    %SWD_output.release.version      = mi_O.SWD_release_version;
    swdOutput.release.version      = swModeOutput.skeleton_version_str;
    swdOutput.release.author       = C.AUTHOR_NAME;
    swdOutput.release.contact      = C.AUTHOR_EMAIL;
    swdOutput.release.institute    = C.INSTITUTE;
    swdOutput.release.modification = swModeOutput.SWD_release_modification;
    swdOutput.release.file         = masterFilename;
    
    % Validate output datasets release version
    % ----------------------------------------
    % RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
    % NOTE: It is hard to thoroughly follow the description, but the end result should be under
    % release_dataset-->version-->pattern (not to be confused with release-->version--pattern).
    if isempty(regexp(swdOutput.release.version, '^[0-9]{2}$', 'once'))
        error('BICAS:get_sw_descriptor:Assertion:IllegalConfiguration', ...
            'Illegal S/W descriptor output release version "%s". This indicates a hard-coded configuration bug.', swdOutput.release.version)
    end
        
    swdMode.outputs.(swModeOutput.JSON_output_file_identifier) = swdOutput;
end



end
