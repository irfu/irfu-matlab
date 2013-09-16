function out=caa_meta(varargin)
% CAA_META return meta data structure
%
%	CAA_META(string1,string2,..) returns name of datasets
%			matching string1,string2,..
%
%	CAA_META(dataset) displays dataset informations
%	out=LOCAL.C_CAA_META(dataset) returns dataset information in structure out
%
%	CAA_META(variableName) displays CAA variableName informations
%	out=LOCAL.C_CAA_META(variableName) returns variable information in structure out
%
%
%	Examples:
%		caa_meta('C1','PEA')
%		caa_meta C1 PEA 
%       d=caa_meta('C4_CP_PEA_MOMENTS')
%       caa_meta('B_Vec_xyz_ISR2__C1_CP_EFW_L2_BB')
%

%	c_caa_meta('create') create file with all structures

%% Defaults and variables
persistent s datasetNames indexFile
indexFileDefault = 'indexCaaMeta_v2'; % default file name, v2 added in 20130618
linkUrlFile = ['http://www.space.irfu.se/cluster/matlab/' ...
	indexFileDefault '.mat'];

%% empty arguments > show help
if nargin==0,
	help caa_meta;
	return;
end

%% Create index file
if nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'create')
	irf_log('fcal','Getting all metadata from CAA, be very patient...');
	urlMetaData = 'http://caa.estec.esa.int/caa_query/?uname=vaivads&pwd=caa&dataset_id=*&metadata=1';
	%urlMetaData = 'http://caa.estec.esa.int/caa_query/?uname=vaivads&pwd=caa&dataset_id=C*_CP_PEA_3DXPH_DEFLUX&metadata=1';
	tempDir = tempname;
	mkdir(tempDir);
	cd(tempDir);
	tempFile = [tempname(tempDir) '.zip'];
	tempFile = urlwrite(urlMetaData,tempFile);
	unzip(tempFile);
	irf_log('fcal','Creating structures');
	d = dir;
	ii=arrayfun(@(x) any(strfind(x.name,'CAA')),d);
	xmlDir = d(ii).name;
	d = dir([xmlDir '/*.XML']);
	isub = [d(:).isdir]; 
	d(isub)=[];
	nameFiles = {d.name}';
	nameDatasetList = cell(numel(nameFiles),1); % needed to construct argument for saving at the end
	datasetList = cell(numel(nameFiles),1);
	for iFile=1:numel(nameFiles)
		fileName=nameFiles{iFile};
		if fileName(1)==' ', disp(['Starts with space:' fileName]);continue;end % in case space in beginning of filename
		nameDataset=fileName(1:find(fileName=='.',1)-1);
		nameDataset(strfind(nameDataset,'-')) = '_';
		datasetList{iFile} = [nameDataset]; % need to add space at end because it is used concatenating for save
		nameDatasetList{iFile} = [nameDataset ' ']; % need to add space at end because it is used concatenating for save
		disp(nameDataset);
		s=xml2struct([xmlDir '/' fileName]);
		if isfield(s,'DATASETS'),
			s=s.DATASETS;
		end
		if isfield(s,'DATASET_METADATA'),
			s=s.DATASET_METADATA;
		end
		meta.(nameDataset) = s; 
	end
	cd ..;
	rmdir(tempDir,'s')
	eval(['save -v7 ' indexFileDefault ' -struct meta ', horzcat(nameDatasetList{:})])
	eval(['save -v7 ' indexFileDefault ' datasetList -append'])
	disp(['Please move ''' indexFileDefault '.mat'' file from the current directory to']);
	disp(linkUrlFile);
	return
end
	
%% Locating index file
if isempty(indexFile) || ~exist(indexFile,'file'), % first usage or file removed
	indexFile = [indexFileDefault '.mat']; % default file name
	if ~exist(indexFile,'file');
		indexFile = irf.get_file(linkUrlFile,'caa',indexFileDefault);
	end
end


if 	nargin==1 && ischar(varargin{1}) && any(strfind(varargin{1},'__')) % CAA variable
	varName=varargin{1};
	dd=regexp(varName, '__', 'split');
	dd{2}(strfind(dd{2},'-')) = '_';	
	dd{2} = upper(dd{2});
	metaData = getfield(load(indexFile,dd{2}),dd{2});
	par=metaData.PARAMETERS.PARAMETER;
	iVar= cellfun(@(x) strcmp(x.PARAMETER_ID.Text,varName),par);
	parVar=par{iVar};
	display_fields(parVar);
	if nargout==1, out=parVar;end
else
	if isempty(s),
		datasetNames = getfield(load(indexFile,'datasetList'),'datasetList');
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
				out = load(indexFile,datasetNames{iSelected});
			elseif any(iEqual)
				out = load(indexFile,datasetNames{iEqual});
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
			dataSet=getfield(load(indexFile,datasetNames{ind}),datasetNames{ind});
			try
				display_fields(dataSet);
				parameters=dataSet.PARAMETERS.PARAMETER;
				disp('PARAMETERS:');
				for j=1:numel(parameters),
					disp([num2str(j) '. ' mat_output(parameters{j}.PARAMETER_ID.Text)]);
				end
			catch
			end
			disp(' ');
		end
		if sum(iSelected) > 1,
			disp('----');
			disp([num2str(sum(iSelected)) ' datasets correspond selection']);
			cellfun(@(x) fprintf('%s\n',x),vertcat(mat_output(datasetNames(iSelected),1)), 'UniformOutput',false);
			return;
		end
		
	end
end
end
%% Functions
function display_fields(dataSet)
disp(' ');
fn=fieldnames(dataSet);
for j=1:numel(fn)
	if isfield(dataSet.(fn{j}),'Text'),
		if strfind(dataSet.(fn{j}).Text,' ') % text includes many words
			disp([fn{j} ': ' dataSet.(fn{j}).Text]);
		else
			disp([fn{j} ': ' mat_output(dataSet.(fn{j}).Text)]);
		end
	end
end
end
function outStr=mat_output(inStr,forceFlag)
% if string is caa variable output link
% if forceFlag defined and equal to one, force matlab linked output
if nargin == 1, 
	forceFlag = false;
end
if ~forceFlag && any(strfind(inStr,'__')),
	forceFlag = true;
end

if forceFlag,
	if ischar(inStr)
		outStr=['<a href="matlab: caa_meta ' ...
			inStr '">' inStr '</a>' ];
	elseif iscell(inStr)
		outStr = cellfun(@(x) ['<a href="matlab: caa_meta ' ...
			x '">' x '</a>' ],inStr, 'UniformOutput',false);
	end
else
	outStr=inStr;
end
end

