function [res,resdataobject,resmat,resunit] = c_caa_var_get(varargin)
%C_CAA_VAR_GET(varName)  load CAA variable varName
%	Dataobject name is derived from the name of varName. If dataobject is in the
%	memory then variable is loaded from there, otherwise it is loaded from the
%	cdf file in CAA directory.
%
%C_CAA_VAR_GET(dobj,varName)  load from dataobject 'dobj'. If dobj is empty 
%	the output is the same as from C_CAA_VAR_GET(varName)
%
%  C_CAA_VAR_GET(varname,'showdep') show dependencies of the variables
%  caa=C_CAA_VAR_GET(varName)         return only in the caa format
%  caa=C_CAA_VAR_GET(varname,'caa')   return only in the caa format
%  var=C_CAA_VAR_GET(varname,'mat')   return only in matlab matrix format
%  dobj=C_CAA_VAR_GET(varname,'dobj') return only data object
%  unit=C_CAA_VAR_GET(varname,'unit') return only units
%
%  [caa,dobj,var,unit]=C_CAA_VAR_GET(varName)
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
if nargin == 0, help c_caa_var_get;return;end
getAllData = true;									% default read all data
getCaa  = true;										% default return variable in caa form
getDobj = false; if nargout>1, getDobj = true;end	% whether to get dataobj
getMat  = false; if nargout>2, getMat  = true;end	% whether to calculate mat variable
getUnit = false; if nargout>3, getUnit = true;end	% whether to get the unit of variable
getMatOnly  = false;
getDobjOnly = false;
getUnitOnly = false;
getFromFile = false;	% reads data from file only if dataobj not in memory
dobjSpecified = false;  % default dataobject is not given as input
dataobject = [];        % default dataobject is empty
args = varargin;        % input argument list
%% Define varName
if isa(args{1},'dataobj'), % dataobject specified as first input argument
	dataobject = args{1};
	dobjSpecified = true;
	args(1)=[];
elseif isempty(args{1}) % dataobject is specified empty
	args(1)=[];         
end

if ~isempty(args)
	varNameList = args{1};
	args(1) = [];
else
	irf.log('critical','not enough input arguments');
	error('not enough input arguments');
end

% Check if ? in varname, then loop through all s/c
if ischar(varNameList) && any(strfind(varNameList,'?'))
	for iSc = 1:4,
		strSc = num2str(iSc);
		tempVarName = strrep(varNameList,'?',strSc);
		res.(['C' strSc]) = c_caa_var_get(dataobject,tempVarName,args{:});
	end
	return
end
while numel(args)
	switch(lower(args{1}))
		case {'show_dependencies','showdep'} % only show dependencies
			if isempty(dataobject)
				dobjName=get_dataobj_name(varNameList);
				existsDobjInCaller = evalin('caller',['exist(''' dobjName ''',''var'')']);
				if existsDobjInCaller,
					dataobject=evalin('caller',dobjName);
				else
					caa_load(dobjName,'nowildcard');
					dataobject=eval(dobjName);
				end
			end
			showdep(dataobject,varNameList)
			return;
		case 'file' % force reading from file
			getFromFile = true;
			args(1)     = [];
		case 'caa' % return caa format only
			% default behaviour, do nothing
			args(1)=[];
		case 'mat' % return matlab format only
			getMatOnly = true;
			getMat     = true;
			getCaa     = true;
			args(1)    = [];
		case 'dobj' % return matlab format only
			getDobjOnly = true;
			getDobj     = true;
			getCaa      = false;
			args(1)     = [];
		case {'unit','units'} % return matlab format only
			getUnitOnly = true;
			getUnit     = true;
			getCaa      = false;
			args(1)     = [];
		case 'tint'                          % load specified time interval
			if numel(args)>1
				if isnumeric(args{2})
					tint = args{2};
				elseif ischar(args{2})
					tint = irf_time(args{2},'iso2tint');
				else
					irf.log('critical','wrongArgType : tint must be numeric or iso')
					error('tint must be numeric or iso');
				end
				args(1:2)=[];
			else
				irf.log('critical','wrongArgType : tint value is missing')
				error('tint value missing');
			end
			getAllData=0;
		otherwise
			irf.log('critical',['Unknown input parameter  ''' args{1} ''' !!!'])
			error('unknown input parameter');
	end
end

%%

if ischar(varNameList),
	varNameList = {varNameList};
end

if ~iscellstr(varNameList)
	irf.log('critical','varName incorrect format');
	error('varName incorrect format');
end

if nargout, % initialize return variables to empty
	res=cell(size(varNameList));
	resdataobject=res;resmat=res;resunit=res;
else
	return;
end

noDataReturned = true;

for iVar = 1: numel(varNameList)
	varName = varNameList{iVar};
	testLocalCaaRepository = false; % default test local CAA directory and not local repository
	if any(strfind(varName,'__'))   % variable name specified as input
		if dobjSpecified
			noDataReturned = false;
		else
			dataobj_name=get_dataobj_name(varName);
			if getAllData && ~getFromFile &&	evalin('caller',['exist(''' dataobj_name ''',''var'')']),
				dataobject=evalin('caller',dataobj_name);
				noDataReturned = false;
				irf.log('warning',[dataobj_name ' exist in memory. NOT LOADING FROM FILE!'])
			else
				if getAllData,
					caa_load(dataobj_name,'nowildcard');
				else
					caa_load(dataobj_name,'tint',tint,'nowildcard');
				end
				if exist(dataobj_name,'var'), % success loading data
					eval(['dataobject=' dataobj_name ';']);
					noDataReturned = false;
				else
					irf.log('warning',[dataobj_name ' could not be loaded!']);
					if (getMat || getCaa || getDobj) && ~getAllData && local.c_read('test')
						testLocalCaaRepository = true;
						irf.log('notice','will test if data are in local CAA data repository.');
					else
						continue;
					end
				end
			end
		end
		if getCaa, % save variable
			if testLocalCaaRepository
				ttt = local.c_read(varName,tint,'caa');
				if isempty(ttt),
					irf.log('warning','NO DATA in repository!');
				end
				noDataReturned = false;
				res{iVar} = ttt;
			else
				res{iVar}=getv(dataobject,varName);
			end
		end
		if getMat % save variable in matlab matrix format
			if testLocalCaaRepository
				ttt = local.c_read(varName,tint,'mat');
				if isempty(ttt),
					irf.log('warning','NO DATA in repository!');
				end
				noDataReturned = false;
				resmat{iVar} = ttt;
			else
				resmat{iVar}=getmat(dataobject,varName);
			end
		end
		if getUnit % save variable unit
			resunit{iVar}=getunits(dataobject,varName);
		end
		if getDobj,% save dataobject
			if testLocalCaaRepository
				ttt = local.c_read(varName,tint,'dobj');
				if isempty(ttt),
					irf.log('warning','NO DATA in repository!');
				end
				noDataReturned = false;
				resdataobject{iVar} = ttt;
			else
				resdataobject{iVar}=dataobject;
			end
		end
	end
end

if noDataReturned, % nothing is loaded, return empty
	irf.log('warning','Nothing is loaded')
	res=[];resdataobject=[];resmat=[];resunit=[];
elseif numel(varNameList) == 1, % return variables and not cell arrays
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


function dataobj_name=get_dataobj_name(varName)
% obtain data object name from full variable name
dd=regexp(varName, '__', 'split');
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
