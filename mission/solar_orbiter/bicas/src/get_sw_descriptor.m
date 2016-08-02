%===================================================================================================
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created ~2016-06-01
%
% Return MATLAB structure "exactly" corresponding to the S/W descriptor specified by the RCS ICD.
% Should effectively return a constant.
%
function SW_descriptor = get_sw_descriptor()
%
% PROPOSAL: Have this function add the prefix "input_" to all modes[].inputs[].input (in SWD)
% and only store the "CLI_parameter_suffix" in the BICAS constants structure instead.
%   CON: The main function also uses the full CLI parameter. It too would have to add the prefix.
%        ==> The same "input_" prefix is specified in two places ==> Bad practice.
%
% PROPOSAL: Implement checks on the information.
%   PROPOSAL: All versions, CLI parameters, dates, dataset IDs on the right format.
%   PROPOSAL: No mode name, dataset ID doubles.
%   PROPOSAL: Check that there are no additional structure fields.
%   PROPOSAL: Check that paths, filenames are valid.
%   NOTE: Already implicitly checks that all the needed fields exist (since they are read here).
%   NOTE: Checks on the main constants structure will (can) only happen if this file is executed, not if 
%         the S/W as a whole is (by default).
%
% PROBLEM: Any validation checks here are not run if bicas_constants.get_constants is called, but
% get_sw_descriptor is not, even if the corresponding information from bicas_constants.get_constants is used.

global CONSTANTS

D.identification = CONSTANTS.C.SWD_identification;
D.release        = CONSTANTS.C.SWD_release;
D.environment    = CONSTANTS.C.SWD_environment;
D.modes = {};

for i = 1:length(CONSTANTS.sw_modes)
    CLI_parameter = CONSTANTS.sw_modes{i}.CLI_parameter;
    
    C_sw_mode = CONSTANTS.get_C_sw_mode_full(CLI_parameter);
    D.modes{end+1} = generate_sw_descriptor_mode(CONSTANTS.C, C_sw_mode);
end

SW_descriptor = D;

end



%===================================================================================================
% Create data structure for a S/W mode corresponding to the information in the JSON S/W descriptor.
%===================================================================================================
function SWD_mode = generate_sw_descriptor_mode(C, C_sw_mode)
%
% Variable naming convention:
%    SWD = S/W descriptor

SWD = [];
SWD.name    = C_sw_mode.CLI_parameter;
SWD.purpose = C_sw_mode.SWD_purpose;

for x = C_sw_mode.inputs
    mi_O = x{1};
    SWD_input = [];
    
    SWD_input.version    = mi_O.dataset_version_str;
    SWD_input.identifier = mi_O.dataset_ID;
    
    SWD.inputs.(mi_O.CLI_parameter) = SWD_input;
end

for x = C_sw_mode.outputs
    mi_O = x{1};
    SWD_output = [];
    
    [master_CDF_path, master_filename] = get_master_CDF_path(mi_O.dataset_ID, mi_O.dataset_version_str);
        
    SWD_output.identifier  = mi_O.dataset_ID;
    SWD_output.name        = mi_O.SWD_name;
    SWD_output.description = mi_O.SWD_description;
    SWD_output.level       = mi_O.SWD_level;
    SWD_output.release.date         = mi_O.SWD_release_date;
    SWD_output.release.version      = mi_O.dataset_version_str;
    SWD_output.release.author       = C.author_name;
    SWD_output.release.contact      = C.author_email;
    SWD_output.release.institute    = C.institute;
    SWD_output.release.modification = mi_O.SWD_release_modification;
    SWD_output.release.file         = master_filename;
        
    SWD.outputs.(mi_O.JSON_output_file_identifier) = SWD_output;
end

SWD_mode = SWD;

end
