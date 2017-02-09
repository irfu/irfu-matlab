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
global CONSTANTS

% SWD = The structure which is used for producing the S/W descriptor (SWD) JSON object string. Its fields (field names)
% should NOT follow variable naming conventions since they influence the JSON object string.
swd = [];
swd.identification = CONSTANTS.C.SWD_IDENTIFICATION;
swd.release        = CONSTANTS.C.SWD_RELEASE;
swd.environment    = CONSTANTS.C.SWD_ENVIRONMENT;
swd.modes = {};

for i = 1:length(CONSTANTS.SW_MODES_INFO_LIST)
    cliParameter = CONSTANTS.SW_MODES_INFO_LIST{i}.CLI_PARAMETER;
    
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



% Create data structure for a S/W mode corresponding to the information in the JSON S/W descriptor.
%
function SwdMode = generate_sw_descriptor_mode(C, ExtendedSwModeInfo)
%
% Variable naming convention:
%    SWD = S/W descriptor

SwdMode = struct;
SwdMode.name    = ExtendedSwModeInfo.CLI_PARAMETER;
SwdMode.purpose = ExtendedSwModeInfo.SWD_PURPOSE;

for OutputInfo = ExtendedSwModeInfo.inputs
    SwModeInputInfo = OutputInfo{1};
    
    SwdInput = struct;
    SwdInput.version    = SwModeInputInfo.SKELETON_VERSION_STR;
    SwdInput.identifier = SwModeInputInfo.DATASET_ID;
    
    SwdMode.inputs.(SwModeInputInfo.CLI_PARAMETER) = SwdInput;
end


for iOutput = 1:length(ExtendedSwModeInfo.outputs)
    OutputInfo    = ExtendedSwModeInfo.outputs{iOutput};

    [~, masterFilename] = bicas.get_master_CDF_path(OutputInfo.DATASET_ID, OutputInfo.SKELETON_VERSION_STR);

    SwdOutputInfo = struct;
    SwdOutputInfo.identifier  = OutputInfo.DATASET_ID;
    SwdOutputInfo.name        = OutputInfo.SWD_NAME;
    SwdOutputInfo.description = OutputInfo.SWD_DESCRIPTION;
    SwdOutputInfo.level       = OutputInfo.SWD_LEVEL;
    SwdOutputInfo.release.date         = OutputInfo.SWD_RELEASE_DATE;
    SwdOutputInfo.release.version      = OutputInfo.SKELETON_VERSION_STR;
    SwdOutputInfo.release.author       = C.AUTHOR_NAME;
    SwdOutputInfo.release.contact      = C.AUTHOR_EMAIL;
    SwdOutputInfo.release.institute    = C.INSTITUTE;
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
