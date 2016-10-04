%===================================================================================================
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created ~2016-06-01
%
% Return MATLAB structure "exactly" corresponding to the S/W descriptor specified by the RCS ICD.
% The values are derived entirely from BICAS constants. The structure is NOT directly incorporated in
% the BICAS constants for extra flexibility.
%
function SW_descriptor = get_sw_descriptor()
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

SWD.identification = CONSTANTS.C.SWD_identification;
SWD.release        = CONSTANTS.C.SWD_release;
SWD.environment    = CONSTANTS.C.SWD_environment;
SWD.modes = {};

for i = 1:length(CONSTANTS.sw_modes)
    CLI_parameter = CONSTANTS.sw_modes{i}.CLI_parameter;    
    
    C_sw_mode = CONSTANTS.get_C_sw_mode_full(CLI_parameter);
    
    SWD.modes{end+1} = generate_sw_descriptor_mode(CONSTANTS.C, C_sw_mode);
end

% Validate S/W release version.
% RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
% NOTE: It is hard to thoroughly follow the description, but the end result should be under
% release-->version-->pattern (not to be confused with release_dataset-->version--pattern).
if isempty(regexp(SWD.release.version, '^(\d+\.)?(\d+\.)?(\d+)$', 'once'))
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal S/W descriptor release version "%s". This indicates a hard-coded configuration bug.', SWD.release.version)
end

SW_descriptor = SWD;

end



%===================================================================================================
% Create data structure for a S/W mode corresponding to the information in the JSON S/W descriptor.
%===================================================================================================
function SWD_mode = generate_sw_descriptor_mode(C, C_sw_mode)
%
% Variable naming convention:
%    SWD = S/W descriptor

global ERROR_CODES

SWD_mode = [];
SWD_mode.name    = C_sw_mode.CLI_parameter;
SWD_mode.purpose = C_sw_mode.SWD_purpose;

for x = C_sw_mode.inputs
    mi_I = x{1};
    SWD_input = [];
    
    SWD_input.version    = mi_I.skeleton_version_str;
    SWD_input.identifier = mi_I.dataset_ID;
    
    SWD_mode.inputs.(mi_I.CLI_parameter) = SWD_input;
end

for x = C_sw_mode.outputs
    mi_O = x{1};
    SWD_output = [];
    
    [master_CDF_path, master_filename] = bicas.get_master_CDF_path(mi_O.dataset_ID, mi_O.skeleton_version_str);
        
    SWD_output.identifier  = mi_O.dataset_ID;
    SWD_output.name        = mi_O.SWD_name;
    SWD_output.description = mi_O.SWD_description;
    SWD_output.level       = mi_O.SWD_level;
    SWD_output.release.date         = mi_O.SWD_release_date;
    %SWD_output.release.version      = mi_O.SWD_release_version;
    SWD_output.release.version      = mi_O.skeleton_version_str;
    SWD_output.release.author       = C.author_name;
    SWD_output.release.contact      = C.author_email;
    SWD_output.release.institute    = C.institute;
    SWD_output.release.modification = mi_O.SWD_release_modification;
    SWD_output.release.file         = master_filename;
    
    % Validate output datasets release version.
    % RCS ICD, iss2rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
    % NOTE: It is hard to thoroughly follow the description, but the end result should be under
    % release_dataset-->version-->pattern (not to be confused with release-->version--pattern).
    if isempty(regexp(SWD_output.release.version, '^[0-9]{2}$', 'once'))
        errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal S/W descriptor output release version "%s". This indicates a hard-coded configuration bug.', SWD_output.release.version)
    end
        
    SWD_mode.outputs.(mi_O.JSON_output_file_identifier) = SWD_output;
end



end
