function caa_export_batch(cl_id,outdir,m_vars,m_vers)
%CAA_EXPORT_BATCH run CAA_EXPORT in a batch script
%
% caa_export_batch(cl_id,sdir,m_vars,m_vers)
%
% M_VARS is a cell array with names of varibles for export
% M_VERS is a corresponding array containing subintervals, file version and
%        data quality information
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if length(m_vars) ~= length(m_vars)
  error('M_VARS and M_VERS must have the same number of elements')
end

sp = pwd;
cd(outdir)

for j=1:length(m_vars)
  v = m_vars{j};
  lev = str2double(v(2));
  caa_vs = v(4:end);

  ints = m_vers{j}{:};
  for ii=1:size(ints,1)
    vers_s = num2str(ints(ii,4));
    if ints(ii,4)<10, vers_s = ['0' vers_s]; end %#ok<AGROW>
    irf_log('save', sprintf('Export : %s -- %s %s V%s Q=%d',...
      epoch2iso(ints(ii,1),1), epoch2iso(ints(ii,1)+ints(ii,2),1), ...
      v, vers_s,ints(ii,3)))
    if caa_export(lev,caa_vs,cl_id,ints(ii,3),vers_s,sp,ints(ii,1),ints(ii,2)) > 0
      error('caa_export returned error')
    end
  end
end

cd(sp)
