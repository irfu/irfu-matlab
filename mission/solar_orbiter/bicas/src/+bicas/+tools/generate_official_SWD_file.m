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
% NOTE: The SWD file should ideally also be validated against the RCS ICD SWD
%       validation schema.
%
%
% ARGUMENTS
% =========
% filePath
%       Optional. If not specified, use <BICAS root dir.>/descriptor.json .
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_official_SWD_file(filePath)
% PROPOSAL: Abolish the --swdescriptor argument.
%   CON: No longer needed.

% Normalize filePath
if nargin == 0
  filePath = bicas.utils.get_SWD_file();
end

Bso = bicas.create_default_BSO();
Bso.make_read_only();

Swml    = bicas.swm.get_SWML(Bso);
JsonSwd = bicas.get_SWD(Swml.List);
strSwd  = bicas.utils.JSON_object_str(JsonSwd);

irf.fs.write_file(filePath, uint8(strSwd(:)))
end
