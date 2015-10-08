function m_update(varargin)
% LOCAL.M_UPDATE update index information in MMS directory
%
%    LOCAL.M_UPDATE update index for all datasets
% 
%    LOCAL.M_UPDATE(datasetName) only update index of dataset "datasetName" 
%
%    LOCAL.M_UPDATE(..,'datadirectory',dataDir) look for data in directory
%    "dataDir". The default data directory is /data/mms unless set by
%    environment variable $DATA_PATH_ROOT.
%
% See also:
%	LOCAL.C_UPDATE (similar function but for Cluster)

%% Defaults
% filterDataSet = false; % default assume datasetName not given as input
%% define local data directory
dataPathRoot = getenv('DATA_PATH_ROOT');
if isempty(dataPathRoot), dataPathRoot='/data/mms';	end; % default MMS data path at IRFU

%% Check inputs
% args = varargin;
% while ~isempty(args)
%   if ischar(args{1}) && strcmpi(args{1},'datadirectory')
%     if numel(args) > 1 && ischar(args{2})
%       dataPathRoot = args{2};
%       args(1:2)=[];
%     else
%       irf.log('warning','data directory not given');
%       args(1)=[];
%     end
%   elseif ischar(args{1}) % given dataset name
%     dataSetFilter = args{1};
%     filterDataSet = true; 
%     args(1)=[];
%   else
%     errStr = 'Unknown input parameter';
%     irf.log('critical',errStr); error(errStr);
%   end
% end
% 
% %% Create dataSetArray to update
% cd(dataPathRoot); %% THONI: WHY CD to this dir?
% tmp = dir(dataPathRoot);
% iDir = [tmp(:).isdir]; % find directories
% dataSetArray = {tmp(iDir).name}';
% dataSetArray(ismember(dataSetArray,{'.','..'})) = []; % remove '.' and '..'
% % dataset name should start with characters mms
% indOkDatasets = cellfun(@(x) (numel(x)>0 && strcmp(x(1:3),'mms')), dataSetArray);
% dataSetArray = dataSetArray(indOkDatasets);
% % If datasetName was given as input remove all other sets.
% if filterDataSet, dataSetArray(~ismember(dataSetArray,{dataSetFilter})) = []; end

% Some files have incorrect FillVal (for instance some hk101 sunpulses) so
% only trust epoch inside of interval 2010-01-01T00:00:00.000000000 to
% 2040-01-01T00:00:00.000000000 (which is well within MMS mission life but
% should discard strange epoch such as 1706-... and other invalid values).
validEpoch = irf.tint('2010-01-01T00:00:00.000000000', ...
  '2040-01-01T00:00:00.000000000');

%% Go through all the datasets to be indexed
for iDataSet = 1:1 %numel(dataSetArray)
  %% list files in data set directory
  dataSetArray={'edp_slow'};
  dataSet = dataSetArray{iDataSet};
  irf.log('warning',['Indexing data set: ' dataSet]);
  %% FIXME: MOVE to proper path...
  cd(['/home/thoni/mms/test_data/mms1/edp/slow']);
  % Locate all mms cdf files for given instrument(-s)
  [unixErr, listFiles] = unix('find . -name ''*mms*.cdf'' | sort');
  if(unixErr)
    errStr = 'Error when trying to list files';
    irf.log('critical', errStr); error(errStr);
  else
    % Spit results into separate cells, one per file.
    listFiles = strsplit(listFiles, '\n');
  end
  if isempty(listFiles)
    irf.log('warning', [dataSet ': no data files']);
    index = [];
  else
    % Pre allocate struct output
    index(1:length(listFiles)) = struct('filename',[],'tstart',[],'tstop',[]);
    ind = 1; % ind used to keep track of which index is written.
    for ii = 1:length(listFiles)
      % One extra "\n" may result in an empty listFiles cell.
      if(isempty(listFiles{ii})), continue; end
      % Get file information from the cdf file
      try
        fileInfo = spdfcdfinfo(listFiles{ii});
      catch
        errStr = ['Cannot get file information from: ', listFiles{ii}];
        irf.log('critical', errStr); warning(errStr); % Should perhaps be error()..
        continue; % Try with next file...
      end
      % Use the fact that the primary epoch variable MUST always be the
      % first variable written to file (files can have many different time
      % series in one single file, but the primary time variable SHOULD
      % always be the first if it it ISTP compliant).
      % KeepEpochAsIs is to ensure it is kept as TT2000 (int64).
      if(~strcmp(fileInfo.Variables{1,4}, 'tt2000'))
        errStr = 'Unexpected first variable, not ISTP compliant cdf file';
        irf.log('critical', errStr); warning(errStr); % Should perhaps be error()...
        continue; % Try with next file
      end;
      try
        epoch = spdfcdfread(listFiles{ii}, 'Variable', fileInfo.Variables{1,1}, 'KeepEpochAsIs', true);
      catch
        errStr = ['Cannot read first variable from file: ', listFiles{ii}];
        irf.log('critical', errStr); warning(errStr); % Should perhaps be error()...
        continue; % Try with next file...
      end
      % All files are not always monotonically increasing in time, so sort
      % the time, run unique to ensure dubblets are removed as well.
      epoch = EpochTT(unique(sort(epoch)));
      % Keep only epoch times within valid epoch interval.
      epoch = epoch(epoch.tlim(validEpoch));
      if(length(epoch)<2)
        % Don't bother with invalid files.
        warnStr = ['TT2000 variable did not contain any valid time interval for file: ', listFiles{ii}];
        irf.log('warning', warnStr); warning(warnStr);
        continue;
      end
      % When we reach this nothing has gone wrong so store file name and
      % interval.
      index(ind).filename = listFiles{ii};
      index(ind).tstart = epoch(1).epoch;
      index(ind).tstop = epoch(end).epoch;
      ind = ind + 1;
    end
    % Remove unused index(), which was preallocated for speed.
    index = index(arrayfun(@(s) ~isempty(s.filename), index));
  end
  %% save index, FIXME: store somewhere good.
  eval(['index_' dataSet '=index;']);
  dirsave(dataSet,['index_' dataSet]);
end

end