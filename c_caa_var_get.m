function [res,resdataobject,resmat,resunit] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(var_name)  load CAA variable from dataobject
%	Dataobject name is derived from the name of the CAA variable.
%
%  var = C_CAA_VAR_GET(varname)
%     get the variable in caa form, dataobject 
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
%  var=C_CAA_VAR_GET(varname,option,'file') force to read dataobject from file
%			even if it is already in memory
%
%  C_CAA_VAR_GET(varname,'tint',tint) get only the interval specified by tint
%  (always reads from file in this case, good option for large files)
%
%   If input is cell array, output is also cell array.
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
%   xm=c_caa_var_get('Differential_Particle_Flux__C3_CP_CIS_HIA_PAD_HS_MAG_IONS_PF','mat');


%% Check input options and set dafaults
getAllData = true;									% default read all data
getCaa  = true;										% default return variable in caa form 
getDobj = false; if nargout>1, getDobj = true;end	% whether to get dataobj
getMat  = false; if nargout>2, getMat  = true;end	% whether to calculate mat variable
getUnit = false; if nargout>3, getUnit = true;end	% whether to get the unit of variable
getMatOnly  = false;
getDobjOnly = false;
getUnitOnly = false;
getFromFile = false;	% reads data from file only if dataobj not in memory
varnames = varargin{1};
args     = varargin(2:end);
nargs    = length(args); % number of additional parameters
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
		case 'file' % force reading from file
			getFromFile = true;			
		case 'caa' % return caa format only
			% default behaviour, do nothing
		case 'mat' % return matlab format only
			getMatOnly=1;
			getMat=1;
			getCaa=0;
		case 'dobj' % return matlab format only
			getDobjOnly=1;
			getDobj=1;
			getCaa=0;
		case {'unit','units'} % return matlab format only
			getUnitOnly=1;
			getUnit=1;
			getCaa=0;
		case 'tint'                          % load specified time interval
			if nargs>1
				if isnumeric(args{2})
					tint = args{2};
					l = 2;
				else irf_log('fcal,','wrongArgType : tint must be numeric')
				end
			else irf_log('fcal,','wrongArgType : tint value is missing')
			end
			getAllData=0;
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
  testLocalCaaRepository = false; % default test local CAA directory and not local repository
  if ischar(var_name) && any(strfind(var_name,'__')) % variable name specified as input
    dataobj_name=get_dataobj_name(var_name);
    if getAllData && ~getFromFile &&	evalin('caller',['exist(''' dataobj_name ''',''var'')']),
      dataobject=evalin('caller',dataobj_name);
      jloaded=jloaded+1;
      irf_log('dsrc',[dataobj_name ' exist in memory. NOT LOADING FROM FILE!'])
    else
      if getAllData,
        caa_load(dataobj_name,'nowildcard');
      else
        caa_load(dataobj_name,'tint',tint,'nowildcard');
      end
      if exist(dataobj_name,'var'), % success loading data
        eval(['dataobject=' dataobj_name ';']);
        jloaded=jloaded+1;
      else
        irf_log('dsrc',[dataobj_name ' could not be loaded!']);
		if (getMat || getCaa || getDobj) && ~getAllData && local.c_read('test')
			testLocalCaaRepository = true;
			irf_log('dsrc','will test if data are in local CAA data repository.');
		else
			continue;
		end
      end
    end
    if getCaa, % save variable
		if testLocalCaaRepository
			ttt = local.c_read(var_name,tint,'caa');
			if isempty(ttt),
				irf_log('dsrc','NO DATA in repository!');
            end
            jloaded = jloaded + 1;
            res{jloaded} = ttt;
		else
			res{jloaded}=getv(dataobject,var_name);
		end
    end
    if getMat % save variable in matlab matrix format
		if testLocalCaaRepository
			ttt = local.c_read(var_name,tint,'mat');
			if isempty(ttt),
				irf_log('dsrc','NO DATA in repository!');
			end
			jloaded = jloaded + 1;
            resmat{jloaded} = ttt;
        else
			resmat{jloaded}=getmat(dataobject,var_name);
		end
    end
    if getUnit % save variable unit
      resunit{jloaded}=getunits(dataobject,var_name);
    end
    if getDobj,% save dataobject 
	  if testLocalCaaRepository
		ttt = local.c_read(var_name,tint,'dobj');
		if isempty(ttt),
		  irf_log('dsrc','NO DATA in repository!');
		end
		jloaded = jloaded + 1;
		resdataobject{jloaded} = ttt;
	  else
      	resdataobject{jloaded}=dataobject; 
	  end
    end
  end
end
if jloaded==0, % nothing is loaded, return empty
  irf_log('load','Nothing is loaded')
  res=[];resdataobject=[];resmat=[];resunit=[];
elseif jloaded ==1, % return variables and not cell arrays
  res=res{1};
  resdataobject=resdataobject{1};
  resmat=resmat{1};
  resunit=resunit{1};
end
if getMatOnly,
  res=resmat; return
end
if getDobjOnly,
  res=resdataobject; return
end
if getUnitOnly,
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

