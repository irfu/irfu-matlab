%
% Automatically generate the official descriptor.json SWD file (as mandated by
% RCS ICD). This function is useful when updating the SWD file before officially
% delivering a new version of BICAS to ROC.
%
% NOTE: Uses default settings which is what BICAS running at ROC should use
% (except for specifying MATLAB executable). Therefore does not need to
% explicitly disable unofficial SWMs (settings SWM.L1-L2_ENABLED,
% SWM.L2-L2_CWF-DSR_ENABLED, SWM.L2-L3_ENABLED),
%
%
% ARGUMENTS
% =========
% filePath
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_official_SWD_file(filePath)
% PROPOSAL: Abolish the --swdescriptor argument.
%   CON: No longer needed.

% Normalize filePath
if nargin == 0
  filePath = fullfile(bicas.utils.get_BICAS_root_path(), 'descriptor.json');
end

Bso = bicas.create_default_BSO();
Bso.make_read_only();

Swml    = bicas.swm.get_SWML(Bso);
JsonSwd = bicas.get_SWD(Swml.List);
strSwd  = bicas.utils.JSON_object_str(JsonSwd, Bso.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));

irf.fs.write_file(filePath, uint8(strSwd(:)))
end
