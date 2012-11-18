function c_update
% LOCAL.C_UPDATE update index information in CAA directory
%
% See also:
%	LOCAL.C_READ

% $Id$
% $Revision$  $Date$


dirCaa='/data/caa/CAA';
cd(dirCaa);
tmp=dir(dirCaa);
iDir = [tmp(:).isdir]; % find directories
dataSetArray = {tmp(iDir).name}';
dataSetArray(ismember(dataSetArray,{'.','..'})) = []; % remove '.' and '..'

for iDataSet=1:numel(dataSetArray)
	%% list files in data set directory
	dataSet=dataSetArray{iDataSet};
	replace_minus_in_cis_names;
	listFiles=dir(dataSet);
	iDir = [listFiles(:).isdir]; %# returns logical vector
	listFiles(iDir)=[];
	%% read in file time intervals
	listFileNames=vertcat(listFiles.name);
	tmp=[listFileNames listFileNames(:,end)]; % add one more column at the end
	tmp(:,end)='=';       % changfe end column to character = (used as separator)
	tmp=tmp';
	tt=textscan(regexprep(tmp(:)','__',' '),'%*[^ ] %4f%2f%2f_%2f%2f%2f_%4f%2f%2f_%2f%2f%2f%*s','delimiter','=');
	%% create index
	index.filename=[repmat([dataSet filesep],numel(listFiles),1) listFileNames];
	index.tstart=irf_time([tt{1} tt{2} tt{3} tt{4} tt{5} tt{6}],'vector2epoch');
	index.tend=irf_time([tt{7} tt{8} tt{9} tt{10} tt{11} tt{12}],'vector2epoch');
	eval(['index_' dataSet '=index;']);
	save('caa',['index_' dataSet],'-append');
end

	function replace_minus_in_cis_names
		if strfind(dataSet,'-'),
			irf_log('dsrc',['Replacing minus signs in dataset: ' dataSet]),
			dataSet=strrep(dataSet,'-','_');
			irf_log('dsrc',['New data set: ' dataSet]);
			listFiles=dir(dataSet);
			iDir = [listFiles(:).isdir]; %# returns logical vector
			listFiles(iDir)=[];
			%% read in file time intervals
			listFileNames=vertcat(listFiles.name);
			if strfind(listFileNames,'-'), % files names to be updated
				cd(dataSet);
				for j=1:size(listFileNames),
					fileName=listFileNames(j,:);
					if strfind(fileName,'-'),
						fileNameNew=strrep(fileName,'-','_');
						evalStr=['!mv ' fileName ' ' fileNameNew];
						irf_log('fcal',evalStr);
						eval(evalStr);
					end
				end
			end
		end
	end
end


