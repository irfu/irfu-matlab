function [res,resdataobject,resmat,resunit] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(var_name)  get CAA variable (if necessary load it)
% 
% if input is cell array, output is also cell array
%
% Usage:
%  var= c_caa_var_get(varargin)
%  [var,dataobject] = c_caa_var_get(varargin)
%  [var,dataobject,variable_matlab_format]=c_caa_var_get(varargin)
%  [var,dataobject,variable_matlab_format,variable_unit]=c_caa_var_get(varargin)
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
%   [~,~,xm]=c_caa_var_get('Differential_Particle_Flux__C3_CP_CIS_HIA_PAD_HS_MAG_IONS_PF');
% give values for each dependency and other information

jloaded=0;
if nargout, % initialize return variables to empty 
    res=[];resdataobject=[];resmat=[];
end

flagmat=0; if nargout>2, flagmat=1;end % whether to calculate mat variable
flagunit=0; if nargout>3, flagunit=1;end % whether to get the unit of variable
for j=1:length(varargin),
  var_name=varargin{j};
  if ischar(var_name) && any(strfind(var_name,'__')) % variable name specified as input
    dd=regexp(var_name, '__', 'split');
    if length(dd)==2, % data object can be properly identifide
        dataobj_name=dd{end};
        if strcmp(dataobj_name,'C3_CP_PEA_'), % the bad case of PEACE
            dataobj_name='C3_CP_PEA_MOMENTS';
        end
    elseif length(dd)==3, % the case of PEACE moments
        if strcmp(dd{3},'MOMENTS'),
            dataobj_name=[dd{2}(1:2) '_CP_PEA_' dd{3}];
        end
    end
    dataobj_name(strfind(dataobj_name,'-'))='_'; % substitute '-' with '_'
    if evalin('caller',['exist(''' dataobj_name ''',''var'')']),
      dataobject=evalin('caller',dataobj_name);
      jloaded=jloaded+1;
      disp('Dataobj exist in memory. NOT LOADING FROM FILE!')
    else
      caa_load(dataobj_name);
      if exist(dataobj_name,'var'), % success loading data
      eval(['dataobject=' dataobj_name ';']);
      jloaded=jloaded+1;
      else
          disp([dataobj_name ' could not be loaded!']);
          continue;
      end
    end
    var=getv(dataobject,var_name);
    if flagmat % construct also matlab format
      varmat=getmat(dataobject,var_name);
    end
    if flagunit % obtain variable unit
      varunit=getunits(dataobject,var_name);
    end
    if jloaded == 1
      res=var;
      resdataobject=dataobject;
      if flagmat, resmat=varmat;end
      if flagunit, resunit=varunit;end
    elseif jloaded == 2
      res={res,var};
      resdataobject={resdataobject,dataobject};
      if flagmat, resmat={resmat,varmat}; end
      if flagunit, resunit={resunit,varunit}; end
    elseif jloaded > 2
      res{jloaded}=var;
      resdataobject{j}=dataobject;
      if flagmat, resmat{j}=varmat;end
      if flagunit, resunit{j}=varunit;end
    end
  end
end
