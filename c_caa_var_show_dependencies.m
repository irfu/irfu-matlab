function c_caa_var_show_dependencies(dobj,var_s)
%SHOWDEP(var_s)  show dependencies for a variable (except time)
%
% SHOWDEP(dobj,var_s)
% SHOWDEP(var_s)
% give values for each dependency and other information

error(nargchk(1,2,nargin))

if nargin == 1
  dd=regexp(dobj, '__', 'split');
  dobjname=dd{end};
  var_s=dobj;clear dobj;
  flag_exist_dobj=evalin('caller',['exist(''' dobjname ''',''var'')']);
  if flag_exist_dobj,
    dobj=evalin('caller',dobjname);
  else
    caa_load(dobjname);
    dobj=eval(dobjname);
  end
end

showdep(dobj,var_s)
