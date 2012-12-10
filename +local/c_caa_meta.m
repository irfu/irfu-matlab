function out=c_caa_meta(varargin)
% LOCAL.C_CAA_META return meta data structure
%
%	LOCAL.C_CAA_META(string1,string2,..) returns name of datasets 
%			matching string1,string2,..
%
%	LOCAL.C_CAA_META(dataset) displays dataset informations 
%	out=LOCAL.C_CAA_META(dataset) returns dataset information in structure out
%
%	LOCAL.C_CAA_META(variableName) displays CAA variableName informations 
%	out=LOCAL.C_CAA_META(variableName) returns variable information in structure out
%
%
%	Examples:
%		LOCAL.C_CAA_META('C1','PEA')
%       d=LOCAL.C_CAA_META('C4_CP_PEA_MOMENTS')
%       LOCAL.C_CAA_META('B_Vec_xyz_ISR2__C1_CP_EFW_L2_BB')
%

% $Id$

%	LOCAL.C_CAA_META('create') create file with all structures

persistent s metaNames datasetNames
indexFile = 'indexCaaMeta.mat';
if isempty(s), % first usage 
	if ~exist(indexFile,'file');
		indexFile = [];
		p=tokenize(path,':');
		for j=1:numel(p),
			if exist([p{j} filesep 'indexCaaMeta.mat'],'file'),
				indexFile = [p{j} filesep 'indexCaaMeta.mat'];
			end
		end
		if isempty(indexFile)
			disp('You do not have indexCaaMat.mat file!')
			disp('Please download it to some directory on your matlab path!');
			disp('<a href="http://www.space.irfu.se/cluster/matlab/indexCaaMeta.mat">http://www.space.irfu.se/cluster/matlab/indexCaaMeta.mat</a>');
			return
		end
	end
end
if nargin==0,
	help local.c_caa_meta;
	return;
end

if nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'create')
	irf_log('fcal','Creating structures');
	d = dir('./*.XML');
	isub = [d(:).isdir]; %# returns logical vector
	d(isub)=[];
	nameFiles = {d.name}';
	delme=[];
	save -v7 index delme;
	for iFile=1:numel(nameFiles)
		fileName=nameFiles{iFile};
		if fileName(1)==' ', disp(['Starts with space:' fileName]);continue;end % in case space in beginning of filename
		nameDataset=fileName(1:find(fileName=='.',1)-1);
		nameDataset(strfind(nameDataset,'-')) = '_';
		disp(nameDataset);
		try
			s=xml2struct(fileName);
		catch
			irf_log('dsrc','Error reading!');
			continue;
		end
		if isfield(s,'DATASETS'),
			s=s.DATASETS;
		end
		if isfield(s,'DATASET_METADATA'),
			s=s.DATASET_METADATA;
		end
		eval(['meta_' nameDataset '=s;']);
		%			save(['meta_' nameDataset],'-v7',['meta_' nameDataset]);
		save('index',['meta_' nameDataset],'-append');
	end
elseif 	nargin==1 && ischar(varargin{1}) && any(strfind(varargin{1},'__')) % CAA variable
	varName=varargin{1};
	dd=regexp(varName, '__', 'split');
	load(indexFile,['meta_' dd{2}]);
	eval(['metaData=meta_' dd{2} ';']);
	par=metaData.PARAMETERS.PARAMETER;
	iVar= cellfun(@(x) strcmp(x.PARAMETER_ID.Text,varName),par);
	parVar=par{iVar};
	display_fields(parVar);
	if nargout==1, out=parVar;end
else
	if isempty(s),
		load(indexFile,'s');
		metaNames={s(:).name};
		datasetNames=arrayfun(@(x) x.name(6:end),s,'uniformoutput',false);
	end
	iSelected = true(numel(datasetNames),1);
	iEqual    = false(numel(datasetNames),1);
	for jInp=1:numel(varargin)
		filter=varargin{jInp};
		if ischar(filter)
			iFind=cellfun(@(x) any(strfind(lower(x),lower(filter))),datasetNames);
			iSelected=iSelected & iFind(:);
			if nargin==1, % check if there is identifcal fit
				iEqual=cellfun(@(x) strcmpi(x,filter),datasetNames);
			end
		end
	end
	if sum(iSelected)>0
		if nargout == 1,
			if sum(iSelected)==1
				load('index',metaNames{iSelected});
				eval(['out = ' metaNames{iSelected} ';']);
			elseif any(iEqual)
				load('index',metaNames{iEqual});
				eval(['out = ' metaNames{iEqual} ';']);
			else
				disp('Output not assigned. There are several choices:');
				disp(vertcat(datasetNames(iSelected)));
				out=[];
			end
			return;
		end
		if sum(iSelected)==1 ||  any(iEqual),
			if any(iEqual),
				ind=find(iEqual);
			else
				ind=find(iSelected);
			end
			s=load(indexFile,metaNames{ind});
			fn=fieldnames(s);
			dataSet=getfield(s,fn{1});
			try
				display_fields(dataSet);
				parameters=dataSet.PARAMETERS.PARAMETER;
				for j=1:numel(parameters),
					disp([num2str(j) '. ' parameters{j}.PARAMETER_ID.Text])
				end
			catch
			end	
			disp(' ');
		end
		if sum(iSelected) > 1,
			disp(vertcat(datasetNames(iSelected)));
			return;
		end

	end
end
function display_fields(dataSet)
fn=fieldnames(dataSet);
for j=1:numel(fn)
	if isfield(dataSet.(fn{j}),'Text'),
		disp([fn{j} ': ' dataSet.(fn{j}).Text]);
	end
end

