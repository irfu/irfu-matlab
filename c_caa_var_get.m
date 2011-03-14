function [res,dataobject] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(var_name)  get variable (except time)
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
% give values for each dependency and other information

for j=1:length(varargin),
  if ischar(varargin{j}) && strfind(varargin{j},'__') % variable name specified as input
    dd=regexp(varargin{j}, '__', 'split');
    dataobj_name=dd{end};
    dataobj_name(strfind(dataobj_name,'-'))='_'; % substitute '-' with '_'
    if evalin('caller',['exist(''' dataobj_name ''',''var'')']),
      dataobject=evalin('caller',dataobj_name);
      disp('Dataobj exist in memory. NOT LOADING FROM FILE!')
    else
      caa_load(dataobj_name);
      eval(['dataobject=' dataobj_name ';']);
    end
    res=getv(dataobject,varargin{j});
  end
end
