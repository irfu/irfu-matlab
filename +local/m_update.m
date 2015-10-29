function m_update(varargin)
% LOCAL.M_UPDATE update index information in MMS directory
%
%    LOCAL.M_UPDATE update index for all datasets for all MMS spacecraft
%
%    LOCAL.M_UPDATE(.., 'scId', {'1', '2', '3'}) only update the specified
%    MMS spacecraft. Any combination ('1', '2', '3' and/or '4') or a single
%    spacecraft.
%
%    LOCAL.M_UPDATE(.., 'datasetName',{'edp', 'fpi'}) only update the
%    specified instrument data sets. In this case "edp" and "fpi".
%
%    LOCAL.M_UPDATE(..,'dataDirectory',dataDir) look for data in directory
%    "dataDir". The default data directory is /data/mms unless set by
%    users environment variable $DATA_PATH_ROOT.
%
% See also:
%	LOCAL.C_UPDATE (similar function but for Cluster)

%% Check inputs
inArg = verify_input(varargin);

% Some files have incorrect FillVal (for instance some hk101 sunpulses) so
% only trust epoch inside of interval 2010-01-01T00:00:00.000000000 to
% 2040-01-01T00:00:00.000000000 (which is well within MMS mission life but
% should discard strange epoch such as 1706-... and other invalid values).
validEpoch = irf.tint('2010-01-01T00:00:00.000000000', ...
  '2040-01-01T00:00:00.000000000');

oldPwd = pwd; % Keep old path.

% Go through each s/c to be indexed
for iSc = 1:numel(inArg.scId)
  % Go through all the datasets to be indexed
  for iDataSet = 1:numel(inArg.datasetName)
    dataSet = inArg.datasetName{iDataSet};
    irf.log('warning',['Indexing data set: ' dataSet ' on MMS',...
      inArg.scId{iSc}]);
    newPath = [inArg.datadirectory, filesep, 'mms', inArg.scId{iSc}, ...
      filesep, dataSet];
    if(~isdir(newPath))
      errStr = ['Not a path: ', newPath];
      irf.log('critical', errStr); error(errStr);
    end
    cd(newPath);
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
      % Remove old cdf files already processed, (ie compare with old index).
      if(exist('irfu_index.mat','file'))
        old_index = load('irfu_index','-mat','index');
        % Some old files found.
        old_index = old_index.index;
        % Keep only files which are still present.
        old_index = old_index(ismember({old_index.filename},listFiles));
        % Identify new files to process.
        listFiles = listFiles(~ismember(listFiles, {old_index.filename}));
        old = true;
      else
        old = false;
      end
      irf.log('warning',['Found ',num2str(length(listFiles)),' new files, will add these to index.']);
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
          irf.log('warning', errStr); warning(errStr); % Should perhaps be error()..
          continue; % Try with next file...
        end
        % Use the fact that the primary epoch variable MUST always be the
        % first variable written to file (files can have many different time
        % series in one single file, but the primary time variable SHOULD
        % always be the first if it it ISTP compliant).
        % KeepEpochAsIs is to ensure it is kept as TT2000 (int64).
        if(~strcmp(fileInfo.Variables{1,4}, 'tt2000'))
          errStr = ['Unexpected first variable, not ISTP compliant cdf file: ' listFiles{ii}];
          irf.log('critical', errStr); warning(errStr); % Should perhaps be error()...
          continue; % Try with next file
        end;
        % Some files have zero records written to epoch. Warn and move on.
        if(fileInfo.Variables{1,3} == 0)
          errStr = ['Empty primary Epoch in cdf file: ', listFiles{ii}];
          irf.log('warning', errStr);
          continue; % Try with next file
        end
        try
          epoch = spdfcdfread(listFiles{ii}, 'Variable', fileInfo.Variables{1,1}, 'KeepEpochAsIs', true);
        catch
          errStr = ['Cannot read first variable from file: ', listFiles{ii}];
          irf.log('critical', errStr); %warning(errStr); % Should perhaps be error()...
          continue; % Try with next file...
        end
        % All files are not always monotonically increasing in time, so sort
        % the time, run unique to ensure dubblets are removed as well.
        epoch = EpochTT(unique(sort(epoch)));
        % Keep only epoch times within valid epoch interval.
        epoch = epoch(epoch.tlim(validEpoch));
        if(length(epoch)<2)
          % Don't bother with invalid files.
          errStr = ['TT2000 variable did not contain any valid time interval for file: ', listFiles{ii}];
          irf.log('critical', errStr); %warning(errStr);
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
    % Combine old and new index list.
    if(old)
      % An old index was read.
      if(~isempty(index))
        % Some new files was found.
        new_index = [old_index, index];
        % Sort the combined list by filename.
        [~, ind_sort] = sort({new_index.filename});
        index = new_index(ind_sort);
      else
        % Only old index, no new files.
        index = old_index;
      end
    end
    % save irfu_index.mat
    save([pwd,filesep,'irfu_index'],'index');
  end
end

% Move back to old path.
cd(oldPwd);

  function inputArg = verify_input(tmpVarargin)
    % Process input arguments
    p = inputParser;
    p.CaseSensitive = false; % match regardless of case. 
    % Default values
    default.scId = {'1','2','3','4'}; % All four MMS spacecrafts
    default.dataSet = {'afg', 'asp1', 'asp2', 'aspoc', 'dfg', 'dsp', ...
      'edi', 'edp', 'epd-eis', 'feeps', 'fields', 'fpi', 'hpca', 'mec', ...
      'scm', 'sdp'}; % All instruments
    dataPathRoot = getenv('DATA_PATH_ROOT');
    if isempty(dataPathRoot)
      dataPathRoot = [filesep,'data',filesep,'mms']; % default MMS path at IRFU
    end;
    default.dataPathRoot = dataPathRoot;
    % Validation functions
    validScId = @(x) assert(all(ismember(x,{'1','2','3','4'})) && ...
      numel(x)==numel(unique(x)) && numel(x)<=4, ...
      'scId should be valid: {"1", "2", "3", "4"}.');
    validDataSet = @(x) assert(all(ismember(x, default.dataSet)) && numel(x)==numel(unique(x)), ...
      'datasetName should be valid MMS instruments.');
    validDatadirectory = @(x) assert(isdir(x), ...
      'DataDirectory should be a directory on your system.');
    % Input arguments, processed in this order of not given as explicit
    % arguments or specified in a struct.
    addOptional(p, 'scId', default.scId, validScId);
    addOptional(p, 'datasetName', default.dataSet, validDataSet);
    addOptional(p, 'datadirectory', default.dataPathRoot, validDatadirectory);
    parse(p, tmpVarargin{:});
    inputArg = p.Results;
  end

end
