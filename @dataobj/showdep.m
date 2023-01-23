function showdep(dobj,var_s)
%SHOWDEP(dobj, var_s)  show dependencies for a variable (except time)
%
% SHOWDEP(dobj, var_s)
% give values for each dependency and other information


% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if ~ischar(var_s), error('VAR_S must be a string'), end

var=getv(dobj,var_s);
disp('_______________');
disp(['VARIABLE : ' var_s]);
disp(var);
for j=0:3 % maximum number of dependencies
  if isfield(var,['DEPEND_' num2str(j)]) % there is dependency
    dep_variable_name = eval(['var.DEPEND_' num2str(j)]);
    disp('_______________');
    disp(['DEPEND_' num2str(j) ': ' dep_variable_name]);
    dep=getv(dobj,dep_variable_name);
    disp(dep);
    disp('Values:');
    disp(num2str([(1:length(dep.data(1,:))); dep.data(1,:)]',3));
  elseif isfield(var,['LABEL_' num2str(j)]) % there is labels
    dep_variable_name = eval(['var.LABEL_' num2str(j)]);
    disp('_______________');
    disp(['LABEL_' num2str(j) ': ' dep_variable_name]);
    dep=getv(dobj,dep_variable_name);
    disp(dep);
    disp('_______________');
    disp('Values:');
    disp(squeeze(dep.data(1,:,:))')
  else
  end
end
