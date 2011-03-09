function [res,dataobj] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(var_name)  get variable (except time)
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
% give values for each dependency and other information

for j=1:length(varargin),
  if ischar(varargin{j}) && strfind(varargin{j},'__') % variable name specified as input
    dd=regexp(varargin{j}, '__', 'split');
    dataobj_name=dd{end};
    if evalin('caller',['exist(''' dataobj_name ''',''var'')']),
      dataobj=evalin('caller',dataobj_name);
      disp('Dataobj exist in memory. NOT LOADING FROM FILE!')
    else
      dataobj=caa_load(dataobj_name);
    end
    res=getv(dataobj,varargin{j});
  end
end
