function [res,resdataobject,resmat,resunit] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(var_name)  get CAA variable (if necessary load it)
%
%  var= C_CAA_VAR_GET(varname)
%     get the variable in caa form
%  [var,dataobject] = C_CAA_VAR_GET(varname)
%     return also dataobject
%  [var,dataobject,variable_matlab_format]=C_CAA_VAR_GET(varname)
%     return variable also in matlab matrix form
%  [var,dataobject,variable_matlab_format,variable_unit]=C_CAA_VAR_GET(varname)
%     return also units
%
%  C_CAA_VAR_GET(varname,'showdep') show dependencies of the variables
%  var=C_CAA_VAR_GET(varname,'mat') return only in matlab matrix format
%  caa=C_CAA_VAR_GET(varname,'caa') equivalent to caa=C_CAA_VAR_GET(varname)
%  dobj=C_CAA_VAR_GET(varname,'dobj') return only data object
%  unit=C_CAA_VAR_GET(varname,'unit') return only units
%
%  C_CAA_VAR_GET(varname,'tint',tint) get only the interval specified by tint
%  (always reads from file in this case, good option for large files)
%
%   If input is cell array, output is also cell array.
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
%   [~,~,xm]=c_caa_var_get('Differential_Particle_Flux__C3_CP_CIS_HIA_PAD_HS_MAG_IONS_PF');

% TODO: add options:
% 'source'(parameters 'file','fast')

%% Check input options
flag_read_all_data=1;                    % default read all data
flagvar=1;                               % default return variable in caa form 
flagdobj=0; if nargout>1, flagdobj=1;end % whether to get dataobj
flagmat=0;  if nargout>2, flagmat=1; end % whether to calculate mat variable
flagunit=0; if nargout>3, flagunit=1;end % whether to get the unit of variable
flag_return_mat_only=0;           % default
flag_return_dobj_only=0;          % default
flag_return_unit_only=0;          % default
varnames=varargin{1};
args=varargin(2:end);
nargs=length(args); % number of additional parameters
while nargs
  l = 1;
  switch(lower(args{1}))
    case {'show_dependencies','showdep'} % only show dependencies
      if nargin == 2
        dobjname=get_dataobj_name(varnames);
        flag_exist_dobj=evalin('caller',['exist(''' dobjname ''',''var'')']);
        if flag_exist_dobj,
          dobj=evalin('caller',dobjname);
        else
          caa_load(dobjname,'nowildcard');
          dobj=eval(dobjname);
        end
      else
        dobj=args{2};
      end
      showdep(dobj,varnames)
      return;
    case 'caa' % return caa format only
        % default behaviour, do nothing
    case 'mat' % return matlab format only
      flag_return_mat_only=1;
      flagmat=1;
      flagvar=0;
    case 'dobj' % return matlab format only
      flag_return_dobj_only=1;
      flagdobj=1;
      flagvar=0;
    case {'unit','units'} % return matlab format only
      flag_return_unit_only=1;
      flagunit=1;
      flagvar=0;
    case 'tint'                          % load specified time interval
      if nargs>1
        if isnumeric(args{2})
          tint = args{2};
          l = 2;
        else irf_log('fcal,','wrongArgType : tint must be numeric')
        end
      else irf_log('fcal,','wrongArgType : tint value is missing')
      end
      flag_read_all_data=0;
    otherwise
      irf_log('fcal',['Unknown input parameter  ''' args{1} ''' !!!'])
      l=1;
  end
  args = args(l+1:end);
  nargs=length(args);
  if isempty(args), break, end
end

%%

jloaded=0;
if nargout, % initialize return variables to empty
  res=cell(1,length(varargin(1)));
  resdataobject=res;resmat=res;resunit=res;
else
  return;  
end

for j=1:length(varargin),
  var_name=varargin{j};
  if ischar(var_name) && any(strfind(var_name,'__')) % variable name specified as input
    dataobj_name=get_dataobj_name(var_name);
    if flag_read_all_data && evalin('caller',['exist(''' dataobj_name ''',''var'')']),
      dataobject=evalin('caller',dataobj_name);
      jloaded=jloaded+1;
      disp('Dataobj exist in memory. NOT LOADING FROM FILE!')
    else
      if flag_read_all_data,
        caa_load(dataobj_name,'nowildcard');
      else
        caa_load(dataobj_name,'tint',tint,'nowildcard');
      end
      if exist(dataobj_name,'var'), % success loading data
        eval(['dataobject=' dataobj_name ';']);
        jloaded=jloaded+1;
      else
        irf_log('dsrc',[dataobj_name ' could not be loaded!']);
        continue;
      end
    end
    if flagvar, % save variable
      res{jloaded}=getv(dataobject,var_name);
    end
    if flagmat % save variable in matlab matrix format
      resmat{jloaded}=getmat(dataobject,var_name);
    end
    if flagunit % save variable unit
      resunit{jloaded}=getunits(dataobject,var_name);
    end
    if flagdobj,% save dataobject 
      resdataobject{jloaded}=dataobject; 
    end
  end
end
if jloaded==0, % nothing is loaded, return empty
  res=[];resdataobject=[];resmat=[];resunit=[];
elseif jloaded ==1, % return variables and not cell arrays
  res=res{1};
  resdataobject=resdataobject{1};
  resmat=resmat{1};
  resunit=resunit{1};
end
if flag_return_mat_only,
  res=resmat; return
end
if flag_return_dobj_only,
  res=resdataobject; return
end
if flag_return_unit_only,
  res=resunit; return
end

function dataobj_name=get_dataobj_name(var_name)
% obtain data object name from full variable name
dd=regexp(var_name, '__', 'split');
if length(dd)==2, % data object can be properly identifide
  dataobj_name=dd{end};
  if strcmp(dataobj_name,'C3_CP_PEA_'), % the bad case of PEACE
    dataobj_name='C3_CP_PEA_MOMENTS';
  elseif strcmp(dataobj_name,'C2_CP_PEA_'), % the bad case of PEACE
    dataobj_name='C2_CP_PEA_MOMENTS';
  elseif strcmp(dataobj_name,'C1_CP_PEA_'), % the bad case of PEACE
    dataobj_name='C1_CP_PEA_MOMENTS';
  elseif strcmp(dataobj_name,'C4_CP_PEA_'), % the bad case of PEACE
    dataobj_name='C4_CP_PEA_MOMENTS';
  end
elseif length(dd)==3, % the case of PEACE moments
  if strcmp(dd{3},'MOMENTS'),
    dataobj_name=[dd{2}(1:2) '_CP_PEA_' dd{3}];
  end
end
dataobj_name(strfind(dataobj_name,'-'))='_'; % substitute '-' with '_'

