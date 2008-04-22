function res = getmat(dobj,var_s)
%GETMAT(dobj, var_s)  get a variable in the Cluster AV format
%
% $Revision$  $Date$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data = getv(dobj,var_s);
dep = getdep(dobj,var_s);

plot_data = double(data.data)';

if isfield(dep,'DEPEND_O'), res = [dep.DEPEND_O plot_data];
else res = plot_data;
end