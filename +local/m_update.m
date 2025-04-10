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
    clear index listFiles old_index
    dataSet = inArg.datasetName{iDataSet};
    irf.log('warning',['Indexing data set: ' dataSet ' on MMS',...
      inArg.scId{iSc}]);
    if(ismember(dataSet,{'defatt','defeph'}))
      ancillary = true;
      % ANCILLARY (ASCII) data
      newPath = [inArg.datadirectory, filesep, 'ancillary', filesep, 'mms',...
        inArg.scId{iSc}, filesep, dataSet];
    else
      ancillary = false;
      % SCIENCE (CDF) data
      newPath = [inArg.datadirectory, filesep, 'mms', inArg.scId{iSc}, ...
        filesep, dataSet];
    end
    if(~isdir(newPath))
      errStr = ['Not a path: ', newPath];
      irf.log('critical', errStr); error(errStr);
    end
    cd(newPath);
    % Locate all MMS files for given dataSet
    if(ancillary)
      [unixErr, listFiles] = unix('find . -name ''*MMS*'' | sort');
    else
      [unixErr, listFiles] = unix('find . -name ''*mms*.cdf'' | sort');
    end
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
      % Remove old files already processed, (ie compare with old index).
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
        % Try to read start from file..
        if(ancillary)
          if(strcmp(dataSet,'defatt'))
            cmd = sprintf('grep COMMENT -A1 %s | tail -n1 | awk ''{print $1}''',...
              listFiles{ii});
          elseif(strcmp(dataSet,'defeph'))
            cmd = sprintf('grep -i km/sec -A1 %s | tail -n1 | awk ''{print $1}''',...
              listFiles{ii});
          end
          [unixErr, startTime] = unix(cmd);
          if(unixErr || isempty(startTime))
            errStr = ['Error when trying to get start time from file: ', listFiles{ii}];
            irf.log('critical', errStr); continue;
          else
            % Try to read stop time
            cmd = sprintf('tail -n2 %s | head -n1 | awk ''{print $1}''',...
              listFiles{ii});
            [unixErr, stopTime] = unix(cmd);
            if(unixErr || isempty(stopTime))
              errStr = ['Error when trying to get stop time from file: ', listFiles{ii}];
              irf.log('critical', errStr); continue;
            else
              % Convert startTime and stopTime to tt2000 int64.
              sss = [irf_time([str2double(startTime(1:4)), str2double(startTime(6:8)); ...
                str2double(stopTime(1:4)), str2double(stopTime(6:8))], ...
                'doy>utc_yyyy-mm-dd'), ['T'; 'T'] ];
              if length(startTime) == 19
                % predatt (time string end with "hh:mm:ss\n") add remaining .mmmuuunnnZ
                sss = [sss, [startTime(10:17); stopTime(10:17)], ['.000000000Z';'.000000000Z']]; %#ok<AGROW>
              else
                % defatt, depeph etc (time string end with "hh:mm:ss.mmm\n" add remaining uuunnnZ
                sss = [sss, [startTime(10:21); stopTime(10:21)], ['000000Z';'000000Z']]; %#ok<AGROW>
              end
              epoch = EpochTT(sss);
              % When this point is reached without error, store result.
              index(ind).filename = listFiles{ii};
              index(ind).tstart = epoch.start.ttns;
              index(ind).tstop = epoch.stop.ttns;
              ind = ind + 1;
            end
          end
        else
          listFiles = remove_corrupted_cdf(listFiles);
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
          EpochId = 1; % Assume it is first variable.
          if(~strcmp(fileInfo.Variables{EpochId,4}, 'tt2000'))
            errStr = ['Not ISTP compliant cdf file: ' listFiles{ii}, '. Trying to locate main Epoch.'];
            irf.log('critical', errStr); warning(errStr); % Should perhaps be error()...
            %continue; % Try with next file
            EpochId = strcmp(fileInfo.Variables(:,1),'Epoch');
            EpochId = find(EpochId,1,'first'); % First Epoch match, if any..
            if(isempty(EpochId))
              errStr = 'No Epoch was identified. Skipping this file.';
              irf.log('critical', errStr); warning(errStr);
              continue
            end
          end
          % Some files have zero records written to epoch. Warn and move on.
          if(fileInfo.Variables{EpochId,3} == 0)
            errStr = ['Empty primary Epoch in cdf file: ', listFiles{ii}];
            irf.log('warning', errStr);
            continue; % Try with next file
          end
          try
            epoch = spdfcdfread(listFiles{ii}, 'Variable', fileInfo.Variables{EpochId,1}, 'KeepEpochAsIs', true);
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
          if(isempty(epoch))
            % Don't bother with invalid files.
            errStr = ['TT2000 variable did not contain any valid time interval for file: ', listFiles{ii}];
            irf.log('critical', errStr); %warning(errStr);
            continue;
          elseif(isscalar(epoch))
            errStr = ['TT2000 variable only contained one valid timestamp for file: ',listFiles{ii}];
            irf.log('warning',errStr);
            tmp_epoch = epoch.epoch;
            clear epoch
            epoch.start.epoch = tmp_epoch;
            epoch.stop.epoch = tmp_epoch;
          end
          % When we reach this nothing has gone wrong so store file name and
          % interval.
          index(ind).filename = listFiles{ii};
          index(ind).tstart = epoch.start.epoch;
          index(ind).tstop = epoch.stop.epoch;
          ind = ind + 1;
        end % if ancillary
      end % for ii = 1:length(listFiles)
    end % if isempty(listFiles)
    % Remove unused index(), which was preallocated for speed.
    index = index(arrayfun(@(s) ~isempty(s.filename), index));
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
  end % for iDataSet
end % for iSc

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
      'scm', 'defatt', 'defeph'}; % All instruments as well as defatt & defeph
    dataPathRoot = getenv('DATA_PATH_ROOT');
    if isempty(dataPathRoot)
      dataPathRoot = [filesep,'data',filesep,'mms']; % default MMS path at IRFU
    end
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

  function listFiles = remove_corrupted_cdf(listFiles)
    % Do not try to read files which are corrupted as this result in
    % problem reading the next files to follow...
    %in bash: cdfdump /data/mms/mms1/asp1/srvy/sitl/beam/2015/08/mms1_asp1_srvy_sitl_beam_20150831_v1.0.1.cdf | less
    % Dumping cdf from "/data/mms/mms1/asp1/srvy/sitl/beam/2015/08/mms1_asp1_srvy_sitl_beam_20150831_v1.0.1.cdf"...
    % Program failed for file: /data/mms/mms1/asp1/srvy/sitl/beam/2015/08/mms1_asp1_srvy_sitl_beam_20150831_v1.0.1.cdf at 1.0...
    % CORRUPTED_V3_CDF: Version 3 CDF is corrupted. (/data/mms/mms1/asp1/srvy/sitl/beam/2015/08/mms1_asp1_srvy_sitl_beam_20150831_v1.0.1.cdf)
    corrupted = './srvy/sitl/beam/2015/08/mms1_asp1_srvy_sitl_beam_20150831_v1.0.1.cdf';
    listFiles = listFiles(~ismember(listFiles, corrupted));
  end

end
