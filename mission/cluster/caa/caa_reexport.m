function caa_reexport(s,QUALITY)
%CAA_REEXPORT  ReExport CEF file
%
% caa_reexport(caa_fname,[QUALITY])
%
% See also CAA_EXPORT, CAA_CEFNAME2SPECS
%

if nargin<2, QUALITY = 3; end
[data_level,caa_vs,cl_id,DATA_VERSION,sp,st,dt] = caa_cefname2specs(s);

disp(['ReExport ' caa_vs ' L' num2str(data_level) ' C' num2str(cl_id) ...
  ' V' DATA_VERSION ' Q=' num2str(QUALITY) ' : '...
  epoch2iso(st,1) '--' epoch2iso(st+dt,1) ])

if caa_export(data_level,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt) ~= 0
  error('Error exporting')
end
