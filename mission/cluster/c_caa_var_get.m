function [res,resdataobject,resmat,resunit] = c_caa_var_get(varargin)
%C_CAA_VAR_GET  load CAA variables
% out=C_CAA_VAR_GET(varName,parameters)
%	Returns variable varName in the format specified by parameters.
%	Dataobject name from which to read variable is derived from the name of
%	varName. If varName is cell array, output is also cell array. The
%	dataobject itself is attempted to obtain in priority order from the
%	following locations:
%		- memory
%		- cdf file in ./CAA directory
%		- using local.c_read()
%		- streaming from CAA data base (if parameters 'mat' and 'tint' specified)
%
% C_CAA_VAR_GET(Dataobj,varName) load from dataobject Dataobj. If Dataobj
%	is empty the output is the same as from C_CAA_VAR_GET(varName).
%
% Output variable types are defined in IRF.DATATYPES.
%  caa=C_CAA_VAR_GET(varName)         returns VariableStruct (in future TSeries)
%  caa=C_CAA_VAR_GET(varname,'caa')   returns VariableStruct
%  caa=C_CAA_VAR_GET(varname,'ts')    returns TSeries object
%  var=C_CAA_VAR_GET(varname,'mat')   returns variableMat
%  dobj=C_CAA_VAR_GET(varname,'dobj') returns DataObject
%  unit=C_CAA_VAR_GET(varname,'unit') returns only units
%
%  [caa,dobj,var,unit]=C_CAA_VAR_GET(varName)
%
%  var=C_CAA_VAR_GET(varname,option,'file') force to read dataobject from
%	file even if it is already in memory.
%
%  C_CAA_VAR_GET(varname,'tint',tint) get data from the time interval
%  specified by tint. If tint is empty get all the data.
%
%  C_CAA_VAR_GET(varname,'showdep') show dependencies of the variables
%
% Example:
%   temp=C_CAA_VAR_GET('Data__C4_CP_PEA_PITCH_SPIN_PSD');
%   xm=c_caa_var_get('Differential_Particle_Flux__C3_CP_CIS_HIA_PAD_HS_MAG_IONS_PF','mat');
%
% See also: IRF.DATATYPES

%% Check input options and set dafaults
if nargin == 0, help c_caa_var_get;return;end
getAllData = true;									% default read all data
getCaa  = true;										% default return variable in caa form
getTs   = false;
getDobj = false; if nargout>1, getDobj = true;end	% whether to get dataobj
getMat  = false; if nargout>2, getMat  = true;end	% whether to calculate mat variable
getUnit = false; if nargout>3, getUnit = true;end	% whether to get the unit of variable
getMatOnly  = false;
getDobjOnly = false;
getUnitOnly = false;
getTsOnly   = false;
getFromFile = false;	% reads data from file only if dataobj not in memory
dobjSpecified = false;  % default dataobject is not given as input
Dataobject = [];        % default dataobject is empty
args = varargin;        % input argument list
testDataStreaming = false; % default do not get data by data streaming
existDataobject = false;
% returnOutputAsCellArray  - set later in code, if true output should be cell array
%% Define varName
if isa(args{1},'dataobj') % dataobject specified as first input argument
	Dataobject = args{1};
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
	for iSc = 1:4
		strSc = num2str(iSc);
		tempVarName = strrep(varNameList,'?',strSc);
		res.(['C' strSc]) = c_caa_var_get(Dataobject,tempVarName,args{:});
	end
	return
end
while numel(args)
	switch(lower(args{1}))
		case {'show_dependencies','showdep'} % only show dependencies
			if isempty(Dataobject)
				dobjName=caa_get_dataset_name(varNameList,'_');
				existsDobjInCaller = evalin('caller',['exist(''' dobjName ''',''var'')']);
				if existsDobjInCaller
					Dataobject=evalin('caller',dobjName);
				else
					caa_load(dobjName,'nowildcard');
					Dataobject=eval(dobjName);
				end
			end
			showdep(Dataobject,varNameList)
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
			getCaa     = false;
			args(1)    = [];
		case 'ts' % return TimeSeries only
			getTsOnly  = true;
			getTs      = true;
			getCaa     = false;
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
				if isempty(args{2})
					getAllData = true;
				else
					if isa(args{2},'GenericTimeArray') && length(args{2})==2
						tintTemp = args{2}.epochUnix;
						tint = tintTemp(:)';
					elseif isnumeric(args{2})
						tint = args{2};
					elseif ischar(args{2})
						tint = irf_time(args{2},'utc>tint');
					else
						irf.log('critical','wrongArgType : tint must be numeric or iso')
						error('tint must be numeric or iso');
					end
					getAllData = false;
				end
				args(1:2)=[];
			else
				irf.log('critical','wrongArgType : tint value is missing')
				error('tint value missing');
			end
		otherwise
			irf.log('critical',['Unknown input parameter  ''' args{1} ''' !!!'])
			error('unknown input parameter');
	end
end

%%

if ischar(varNameList)
	returnOutputAsCellArray = false;
	varNameList = {varNameList};
elseif iscellstr(varNameList)
	returnOutputAsCellArray = true;
else
	irf.log('critical','varName incorrect format');
	error('varName incorrect format');
end

if nargout % initialize return variables to empty
	res=cell(size(varNameList));
	resdataobject=res;resmat=res;resunit=res;rests=res;
else
	return;
end

isDataReturned = false;

datasetNameList = caa_get_dataset_name(varNameList,'_');
[datasetNameUniqueList,~,indDatasetUniqueToName] ...
	= unique(datasetNameList);
for iDataset = 1:numel(datasetNameUniqueList)
	dataobjName = datasetNameUniqueList{iDataset};
	indVarNameList = find(indDatasetUniqueToName == iDataset);
	testLocalCaaRepository = false; % default test local CAA directory and not local repository
	if dobjSpecified
		isDataReturned = true;
		existDataobject = true;
	else
		if getAllData && ~getFromFile ...
				&& evalin('caller',['exist(''' dataobjName ''',''var'')'])...
				&& evalin('caller',['isa(' dataobjName ',''dataobj'')'])
			Dataobject=evalin('caller',dataobjName);
			isDataReturned = true;
			existDataobject = true;
			irf.log('warning',[dataobjName ' exist in memory. NOT LOADING FROM FILE!'])
		else
			if getAllData
				caa_load(dataobjName,'nowildcard');
			else
				caa_load(dataobjName,'tint',tint,'nowildcard');
			end
			if exist(dataobjName,'var') % success loading data
				eval(['Dataobject=' dataobjName ';']);
				isDataReturned = true;
				existDataobject = true;
			else
				irf.log('warning',[dataobjName ' could not be loaded!']);
				if (getTs || getMat || getCaa || getDobj) && ~getAllData && local.c_read('test')
					testLocalCaaRepository = true;
					irf.log('notice','will test if data are in local CAA data repository.');
				end
				if getMat && ~getAllData
					if datastore('irfu_matlab','okCeflib')
						testDataStreaming = true;
					end
				end
			end
		end
	end
	varTmpList = varNameList(indVarNameList);
	if getCaa % save variable
		if existDataobject
			res(indVarNameList) = cellfun(@(x) getv(Dataobject,x),varTmpList,'uniformoutput',false);
		elseif testLocalCaaRepository
			ttt = local.c_read(varTmpList,tint,'caa');
			if isempty(ttt)
				irf.log('warning','NO DATA in repository!');
			end
			isDataReturned = true;
			res(indVarNameList) = ttt;
		end
	end
	if getMat % save variable in matlab matrix format
		if existDataobject
			resmat(indVarNameList) = cellfun(@(x) getmat(Dataobject,x),varTmpList,'uniformoutput',false);
		else
			if testLocalCaaRepository
				ttt = local.c_read(varTmpList{1},tint,'mat');
				if isempty(ttt)
					irf.log('warning','NO DATA in local repository!');
				else
					resmat{indVarNameList(1)} = ttt;
					for j=2:numel(varTmpList)
						ttt = local.c_read(varTmpList{j},tint,'mat');
						resmat{indVarNameList(j)} = ttt;
					end
					isDataReturned = true;
				end
			end
			if ~isDataReturned && testDataStreaming
				ttt = c_caa_cef_var_get(varTmpList,'tint',tint,'stream');
				if isempty(ttt)
					irf.log('warning','NO DATA in CAA to stream!');
				else
					resmat(indVarNameList) = ttt;
					isDataReturned = true;
				end
			end
		end
	end
	if getTs % get variable in TimeSeries format
		if existDataobject
			rests(indVarNameList) = cellfun(@(x) get_ts(Dataobject,x),varTmpList,'uniformoutput',false);
		else
			if testLocalCaaRepository
				ttt = local.c_read(varTmpList{1},tint,'ts');
				if isempty(ttt)
					irf.log('warning','NO DATA in local repository!');
				else
					rests{indVarNameList(1)} = ttt;
					for j=2:numel(varTmpList)
						ttt = local.c_read(varTmpList{j},tint,'ts');
						rests{indVarNameList(j)} = ttt;
					end
					isDataReturned = true;
				end
			end
			if ~isDataReturned && testDataStreaming
				ttt = c_caa_cef_var_get(varTmpList,'tint',tint,'stream','ts');
				if isempty(ttt)
					irf.log('warning','NO DATA in CAA to stream!');
				else
					rests(indVarNameList) = ttt;
					isDataReturned = true;
				end
			end
		end
	end
	if getUnit % save variable unit TODO: implement local.c_read and streaming
		if existDataobject
			resunit(indVarNameList) = cellfun(@(x) getunits(Dataobject,x),varTmpList,'uniformoutput',false);
		end
	end
	if getDobj% save dataobject TODO: call to stream data object
		for iVar = indVarNameList(:)'
			varName = varNameList{iVar};
			if existDataobject
				resdataobject{iVar}=Dataobject;
			elseif testLocalCaaRepository
				ttt = local.c_read(varName,tint,'dobj');
				if isempty(ttt)
					irf.log('warning','NO DATA in repository!');
				end
				isDataReturned = true;
				resdataobject{iVar} = ttt;
			end
		end
	end
end

if ~isDataReturned % nothing is loaded, return empty
	irf.log('warning','Nothing is loaded')
	res=[];resdataobject=[];resmat=[];resunit=[];
elseif numel(varNameList) == 1 && ~returnOutputAsCellArray % return variables and not cell arrays
	res=res{1};
	resdataobject=resdataobject{1};
	resmat=resmat{1};
	resunit=resunit{1};
	rests=rests{1};
end
if getMatOnly
	res=resmat; return
end
if getTsOnly
	res=rests; return
end
if getDobjOnly
	res=resdataobject; return
end
if getUnitOnly
	res=resunit; return
end
