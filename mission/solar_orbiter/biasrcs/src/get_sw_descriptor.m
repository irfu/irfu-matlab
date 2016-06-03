% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created ~2016-06-01
%
% Get MATLAB structure corresponding to the SW descriptor specified by the RCS ICD.
%
function SW_descriptor = get_sw_descriptor()

persistent D
if ~isempty(D)
    SW_descriptor = D;
    return
end



C = biasrcs_constants;

D.identification = C.SWD_identification;
D.release        = C.SWD_release;
D.environment    = C.SWD_environment;
D.modes = {};

for mi = C.sw_modes
    D.modes{end+1} = generate_sw_descriptor_mode(C, mi{1});
end

SW_descriptor = D;

end

%===================================================================================================

% mi = mode info
% sw_descriptor_mode = Mode information that can be interpreted as a mode in the SW descriptor (JSON).
function SWD_mode = generate_sw_descriptor_mode(C, mi)
% Variable naming convention:
%    swd = SW descriptor
%    mi = mode info (information stemming from the function argument)

swd = [];
swd.name    = mi.CLI_parameter;
swd.purpose = mi.SWD_purpose;

for x = mi.inputs
    mi_O = x{1};
    swd_input = [];
    
    swd_input.version    = mi_O.dataset_version_str;
    swd_input.identifier = mi_O.dataset_ID;
    
    swd.inputs.(mi_O.CLI_parameter_name) = swd_input;
end

for x = mi.outputs
    mi_O = x{1};
    swd_output = [];
        
    swd_output.identifier  = mi_O.dataset_ID;
    swd_output.name        = mi_O.SWD_name;
    swd_output.description = mi_O.SWD_description;
    swd_output.level       = mi_O.SWD_level;
    swd_output.author      = C.author_name;
    swd_output.release.date         = mi_O.SWD_release_date;
    swd_output.release.version      = mi_O.dataset_version_str;
    swd_output.release.contact      = C.author_email;
    swd_output.release.institute    = C.institute;
    swd_output.release.modification = mi_O.SWD_release_modification;
    swd_output.release.file         = mi_O.master_cdf_filename;
        
    swd.outputs.(mi_O.JSON_output_file_identifier) = swd_output;
end

SWD_mode = swd;

end
