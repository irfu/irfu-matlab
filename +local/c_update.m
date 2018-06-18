function c_update(varargin)
% LOCAL.C_UPDATE update index information in CAA directory
%
%    LOCAL.C_UPDATE update index for all datasets
% 
%    LOCAL.C_UPDATE(datasetName) update only the index datasetName or all datasets beginning with dataSetName. (Ex 'C2_')
%
%    LOCAL.C_UPDATE(..,'datadirectory',dataDir) look for data in directory
%    dataDir. The default data directory is /data/caalocal unless set by
%		datastore('caa','localDataDirectory','/new/path')
%
% See also:
%	LOCAL.C_READ

%% Defaults
filterDataSet = false; % default assume datasetName not given as input
%% define local data directory
dirCaa = datastore('caa','localDataDirectory');
if isempty(dirCaa)
	dirCaa='/data/caalocal';	% default
end
%% Check inputs
args=varargin;
while ~isempty(args)
	if ischar(args{1}) && strcmpi(args{1},'datadirectory')
		if numel(args) > 1 && ischar(args{2})
			dirCaa = args{2};
			args(1:2)=[];
		else
			irf.log('warning','data directory not given');
			args(1)=[];
		end
	elseif ischar(args{1}) % given dataset name
		dataSetFilter = args{1};
		filterDataSet = true; 
		args(1)=[];
	else
		errStr = 'local.c_update: unknown input parameter';
		irf.log('critical',errStr);
		disp(args{1});
		error('local:c_update:input',errStr);
	end
end
%% Create dataSetArray to update
cd(dirCaa);
tmp=dir(dirCaa);
iDir = [tmp(:).isdir]; % find directories
dataSetArray = {tmp(iDir).name}';
dataSetArray(ismember(dataSetArray,{'.','..'})) = []; % remove '.' and '..'
% dataset name should start with character C
indOkDatasets = cellfun(@(x) (numel(x)>0 && (x(1) == 'C')), dataSetArray);
dataSetArray = dataSetArray(indOkDatasets);
if filterDataSet
	dataSetArray=dataSetArray(strmatch({dataSetFilter},dataSetArray));
%	dataSetArray(~ismember(dataSetArray,{dataSetFilter})) = [];
end

%% go through all the datasets to be indexed
for iDataSet=1:numel(dataSetArray)
	%% list files in data set directory
	dataSet=dataSetArray{iDataSet};
	irf.log('warning',['Indexing data set: ' dataSet]);
	dataSet=replace_minus_in_cis_names(dataSet);
	listFiles=dir([dataSet '/*.cdf']);
	iDir = [listFiles(:).isdir]; %# returns logical vector
	listFiles(iDir)=[];
	if numel(listFiles)==0 
		irf.log('warning',[dataSet ': no data files']);
		index=[];
	else
		%% read in file time intervals
		listFileNames=vertcat(listFiles.name);
		tmp=[listFileNames listFileNames(:,end)]; % add one more column at the end
		tmp(:,end)='=';       % change end column to character = (used as separator)
		tmp=tmp';
		tt=textscan(regexprep(tmp(:)','__',' '),'%*[^ ] %4f%2f%2f_%2f%2f%2f_%4f%2f%2f_%2f%2f%2f_V%6s%*s','delimiter','=');
		%% create index
		index.filename=[repmat([dataSet filesep],numel(listFiles),1) listFileNames];
		index.tstart=irf_time([tt{1} tt{2} tt{3} tt{4} tt{5} tt{6}],'vector>epoch');
		index.tend=irf_time([tt{7} tt{8} tt{9} tt{10} tt{11} tt{12}],'vector>epoch');
		index.versionFile=tt{13};
	end
	%% save index
	eval(['index_' dataSet '=index;']);
	dirsave(dataSet,['index_' dataSet]);
end

	function dataSetOut=replace_minus_in_cis_names(dataSet)
		if strfind(dataSet,'CIS') %#ok<STRIFCND>
			if strfind(dataSet,'-') %#ok<STRIFCND>
				irf.log('debug',['Replacing minus signs in dataset: ' dataSet]),
				dataSetNew=strrep(dataSet,'-','_');
				irf.log('debug',['New data set: ' dataSetNew]);
				if ~exist(dataSetNew, 'dir')
					movefile(dataSet,dataSetNew);
				else
					movefile([dataSet '/*.cdf'],dataSetNew);
					rmdir(dataSet); %% can fail that is ok if dir is not empty
				end
				dataSet=dataSetNew;
			end
			listFiles=dir([dataSet '/*.cdf']);
			iDir = [listFiles(:).isdir]; %# returns logical vector
			listFiles(iDir)=[];
			%% read in file time intervals
			listFileNames=vertcat(listFiles.name);
			if strfind(listFileNames(:)','-')  %#ok<STRIFCND> % files names to be updated
				cd(dataSet);
				for j=1:size(listFileNames)
					fileName=listFileNames(j,:);
					if strfind(fileName,'-') %#ok<STRIFCND>
						fileNameNew=strrep(fileName,'-','_');
						irf.log('debug',['moving: ' fileName ' >> ' fileNameNew]);
						movefile(fileName,fileNameNew);
					end
				end
				cd('..');
			end
		end
		dataSetOut=dataSet;
