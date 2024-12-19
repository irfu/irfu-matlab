function out=caa_download(varargin)
% LOCAL.CAA_DOWNLOAD download full datasets from CAA
% In case data already exists on disk, downloads only newer version files.
% Dataset location - default is /data/caalocal but can be specified as
% input parameter or set by datastore('caa','localDataDirectory','/new/directory/location')
% Dataset request timetable information - directory matCaaRequests
% Index location - inside the directory of each dataset
% Current request time table is accessible in variable TTRequest
%
%   LOCAL.CAA_DOWNLOAD(dataset) download all dataset
%
%   LOCAL.CAA_DOWNLOAD(dataset,'indexStart',indexStart) start from the
%   specified index indexStart
%
%   LOCAL.CAA_DOWNLOAD(dataset,'indexList',indexList) download only indexes
%   specified in indexList
%
%   LOCAL.CAA_DOWNLOAD(dataset,'stream') stream the data
%
%   LOCAL.CAA_DOWNLOAD(dataset,'daily') download daily files
%   LOCAL.CAA_DOWNLOAD(dataset,'monthly') download monthly files
%
%   LOCAL.CAA_DOWNLOAD(dataset,'DataDirectory',dataDir) use dataDir as
%   location for data (default dataDir is '/data/caalocal')
%
%   LOCAL.CAA_DOWNLOAD(dataset,'simulate') show statistics on how many
%   intervals would be deleted and downloaded but do not do anything.
%
%   LOCAL.CAA_DOWNLOAD(dataset,'email') send email when finnished
%
%   LOCAL.CAA_DOWNLOAD(...,inputParamCaaDownload) any unrecognized input
%   parameter is parsed to caa_download when downloading data. See help
%   caa_download.
%
%   TTrequest = LOCAL.CAA_DOWNLOAD(...) return resulting time table TTrequest, where
%
%		TTRequest.UserData.Status = 1 - downloaded, 0 - submitted, empty - not processed
%		TTRequest.UserData.Downloadfile zip file to download (important if status=0)
%		TTRequest.UserData.TimeOfRequest
%		TTRequest.UserData.TimeOfDownload
%		TTRequest.UserData.NumberOfAttemptsToDownload
%		TTRequest.UserData.dataset
%		TTRequest.UserData.number - number of entries
%		TTRequest.UserData.version - version of dataset
%
%   LOCAL.CAA_DOWNLOAD(TTrequest) process request time table TTRequest
%
% Dataset downloading takes long time. You will be informed by email when it is
% finished if you have specified your information in datastore. Just execute in
% matlab (substitute Your Name to your name):
%		datastore('local','name','Your Name')
%
% Example:
%		local.caa_download('C1_CP_PEA_MOMENTS')
%
% 	See also CAA_DOWNLOAD

%% Defaults
dataDir = datastore('caa','localDataDirectory');
if isempty(dataDir)
  dataDir			= '/data/caalocal';
end
maxSubmittedJobs      = 10;
maxNumberOfAttempts   = 30;
isInputDatasetName    = false;
sendEmailWhenFinished = false;
streamData            = false; % download cdf files asynchronously
indexStart            = 1;
inputParamCaaDownload = {};
doSimulateDownload    = false; % takes care of flag 'simulate'
doDailyFileDownload   = false; % default is to go by inventory time
doMonthlyFileDownload = false; % default is to go by inventory time

%% check input: get inventory and construct time table if dataset
if nargin == 0
  help local.caa_download
  return
end
args=varargin;
if ischar(varargin{1})
  dataSet=varargin{1};
  if isempty(dataSet)
    irf.log('warning','local.caa_download: dataset name empty!')
    return;
  end
  isInputDatasetName = true;
elseif isa(varargin{1},'irf.TimeTable')
  TTRequest=varargin{1};
  assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
  dataSet=dataset_mat_name(TTRequest.UserData(1).dataset);
else
  irf.log('critical','See syntax: help local.c_caa_download');
  return;
end
args(1)=[];

while ~isempty(args)
  if ischar(args{1}) && strcmpi(args{1},'stream')
    streamData = true;
    irf.log('notice','Streaming data.');
    args(1) = [];
  elseif ischar(args{1}) && strcmpi(args{1},'simulate')
    irf.log('notice','Only simulate the download');
    doSimulateDownload = true;
    args(1) = [];
  elseif ischar(args{1}) && strcmpi(args{1},'email')
    sendEmailWhenFinished = true;
    args(1) = [];
  elseif ischar(args{1}) && strcmpi(args{1},'daily')
    irf.log('notice','Download as daily files');
    doDailyFileDownload = true;
    args(1) = [];
  elseif ischar(args{1}) && strcmpi(args{1},'monthly')
    irf.log('notice','Download as monthly files');
    doMonthlyFileDownload = true;
    args(1) = [];
  elseif ischar(args{1}) && strcmpi(args{1},'indexstart')
    if numel(args) > 1 && isnumeric(args{2})
      indexStart = args{2};
      args(1:2)=[];
    else
      irf.log('warning','indexstart not given');
      args(1)=[];
    end
  elseif ischar(args{1}) && strcmpi(args{1},'indexlist')
    if numel(args) > 1 && isnumeric(args{2})
      indexList = args{2};
      args(1:2)=[];
    else
      irf.log('warning','indexList not given');
      args(1)=[];
    end
  elseif ischar(args{1}) && strcmpi(args{1},'datadirectory')
    if numel(args) > 1 && ischar(args{2})
      dataDir = args{2};
      args(1:2)=[];
    else
      irf.log('warning','data directory not given');
      args(1)=[];
    end
  else
    inputParamCaaDownload{end+1} = args{1}; %#ok<AGROW>
    args(1) = [];
  end
end
%% If sending email check that all servers and addresses defined
% datastore('local',..) keeps info on email and servers
if sendEmailWhenFinished
  if exist('sendmail','file')==2
    sendEmailFrom = datastore('local','sendEmailFrom');
    if isempty(sendEmailFrom)
      disp('Please define email from which MATLAB should send email');
      disp('This should be your local email where you are running MATLAB.');
      disp('For example in IRFU: username@irfu.se');
      disp('Execute in matlab (adjust accordingly):')
      disp(' ');
      disp('>  datastore(''local'',''sendEmailFrom'',''username@irfu.se'')');
      disp(' ');
      return;
    end
    sendEmailSmtp = datastore('local','sendEmailSmtp');
    if isempty(sendEmailSmtp)
      disp('Please define your local SMTP server. ');
      disp('For example in IRFU: sol.irfu.se');
      disp('Execute in matlab (adjust accordingly):')
      disp(' ');
      disp('>  datastore(''local'',''sendEmailSmtp'',''sol.irfu.se'')');
      disp(' ');
      return;
    end
    sendEmailTo = datastore('local','email');
    if isempty(sendEmailTo)
      disp('Please specify your email.');
      disp('For example: name@gmail.com');
      disp('Execute in matlab (adjust accordingly):')
      disp(' ');
      disp('>  datastore(''local'',''email'',''name@gmail.com'')');
      disp(' ');
      return;
    end
  else
    disp('Your system does not allow sending emails');
    disp('"sendmail" command does not exist.')
    sendEmailWhenFinished = false;
  end
end
%% change to data directory
if exist(dataDir,'dir')
  disp(['!!!! Changing directory to ' dataDir ' !!!']);
  cd(dataDir);
else
  irf.log('critical',['DOES NOT EXIST data directory: ' dataDir ' !!!']);
  out=[];
  return;
end

%% Check inventory when input is dataset name
%  check also if daily files are downloaded
if isInputDatasetName
  irf.log('warning','Checking list of available times');
  if doMonthlyFileDownload
    TT=caa_download(['list:' dataSet]);
    if numel(TT)==0
      disp('Dataset does not exist or there are no data');
      return;
    end
    tminDatenum = irf_time(TT.TimeInterval(1),'epoch>vector');
    tmaxDatenum = irf_time(TT.TimeInterval(2),'epoch>vector');
    startYear   = tminDatenum(1);
    startMonth  = tminDatenum(2);
    endYear     = tmaxDatenum(1);
    endMonth    = tmaxDatenum(2)+1;
    if endMonth == 13
      endYear  = endYear + 1;
      endMonth = 1;
    end
    tVec = zeros((1 + endYear - startYear)*12,1);
    for iYear = startYear:endYear
      for iMonth = 1:12
        tVec((iYear-startYear)*12+iMonth) = irf_time([iYear iMonth 1 0 0 0],'vector>epoch');
      end
    end
    tStart      = tVec(startMonth  :end-(12-endMonth)-1);
    tEnd        = tVec(startMonth+1:end-(12-endMonth)  );
    TTRequest   = irf.TimeTable([tStart tEnd]);
  elseif doDailyFileDownload
    TT=caa_download(['list:' dataSet]);
    if numel(TT)==0
      disp('Dataset does not exist or there are no data');
      return;
    end
    tminDatenum = irf_time(TT.TimeInterval(1),'epoch>datenum');
    tmaxDatenum = irf_time(TT.TimeInterval(2),'epoch>datenum');
    tminDatenum = floor(tminDatenum);
    tmaxDatenum = floor(tmaxDatenum) + 1;
    tStart      = irf_time((tminDatenum : tmaxDatenum)'  ,'datenum>epoch');
    tEnd        = irf_time((tminDatenum : tmaxDatenum)'+1,'datenum>epoch');
    TTRequest   = irf.TimeTable([tStart tEnd]);
  else
    TT=caa_download(['inventory:' dataSet]);
    if numel(TT)==0
      disp('Dataset does not exist or there are no data');
      return;
    end
    TTRequest=TT;
  end
  assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
end
%% check which time intervals are already downloaded, remove obsolete ones
dataSetDir = [dataDir filesep dataSet];
indNewIntervals = true( numel(TTRequest),1); % default
if exist(dataSetDir,'dir')
  % read index of already present dataset files
  indexDataSetName = ['index_' dataSet];
  indexDataSetFileName = [dataSetDir filesep indexDataSetName  '.mat'];
  if ~exist(indexDataSetFileName,'file')
    local.c_update(dataSet);
  end
  dirload(dataSetDir,indexDataSetName);
  indexDataSet = eval(indexDataSetName);
  if isempty(indexDataSet)
    TTindex = [];
  else
    TTindex = irf.TimeTable([indexDataSet.tstart indexDataSet.tend]);
    lastVersion = 0;
    for ii = numel(indexDataSet.versionFile):-1:1
      version = str2double(indexDataSet.versionFile{ii});
      if version > lastVersion, lastVersion = version; end
      TTindex.UserData(ii).version  = version;
      TTindex.UserData(ii).filename = indexDataSet.filename(ii,:);
    end
  end
  % define request
  if numel(TTindex)>0
    irf.log('warning','Previous data exist, merging...');
    TTindex=sort(TTindex);

    % obtain file list that are ingested since the last data file
    TTfileList=caa_download(['fileinventory:' dataSet]);
    indNewFiles = false(1,numel(TTfileList));
    for j = 1:numel(indNewFiles)
      newVersion = str2double(TTfileList.UserData(j).caaIngestionDate([3 4 6 7 9 10]));
      if newVersion > lastVersion
        indNewFiles(j) = true;
        irf.log('debug',['Ingested since last #' num2str(j) ...
          ' time interval: ' irf_time(TTfileList.TimeInterval(j,:),'tint>utc_yyyy-mm-dd')]);
      end
    end
    TTfileList = select(TTfileList,find(indNewFiles));

    % find TTRequest intervals with new files in them
    [~,ii]=overlap(TTRequest,TTfileList);
    indNewIntervals(:)  = false;
    indNewIntervals(ii) = true;

    % check which old intervals to be removed/updated
    tintInd = TTindex.TimeInterval;
    tintReq = TTRequest.TimeInterval;
    indOldObsoleteIntervals = false(numel(TTindex),1);
    indOldToUpdateIntervals = false(numel(TTindex),1);
    nOldIntervals = numel(TTindex);
    jIndex        = nOldIntervals;
    for iReq=numel(TTRequest):-1:1
      if tintInd(jIndex,2) < tintReq(iReq,1) % new time interval
        continue;
      else
        while(tintInd(jIndex,2)>tintReq(iReq,1))
          if	abs(tintInd(jIndex,2)-tintReq(iReq,2))<=5 && ...% interval comparison to 5s (1spin) precision
              abs(tintInd(jIndex,1)-tintReq(iReq,1))<=5
            if indNewIntervals(iReq) % new files available for interval
              indOldToUpdateIntervals(jIndex) = true;
              irf.log('debug',['Old interval to update #' num2str(jIndex) ' '...
                irf_time(tintInd(jIndex,:),'tint>utc')]);
            end
          else
            indOldObsoleteIntervals(jIndex) = true;
            irf.log('debug',['Obsolete interval #' num2str(jIndex) ' '...
              irf_time(tintInd(jIndex,:),'tint>utc')]);
          end
          jIndex = jIndex-1;
          if jIndex==0, break; end
        end
      end
      if jIndex==0, break; end
    end
    nOldObsoleteIntervals = sum(indOldObsoleteIntervals);
    nOldToUpdateInterval  = sum(indOldToUpdateIntervals);
    % Removing old obsolete files
    if ~doSimulateDownload
      remove_datafiles(TTindex,indOldObsoleteIntervals,dataDir);
      local.c_update(dataSet); % so that index is as soon as files are removed
    end
    % find new intervals that do not overlap with old ones
    TTindexUnchanged = select(TTindex,~indOldObsoleteIntervals & ~indOldToUpdateIntervals);
    [~,ii]=overlap(TTRequest,TTindexUnchanged);
    indNonOverlap = ~(false.*indNewIntervals);
    indNonOverlap(ii) = false;
    indNewIntervals(indNonOverlap) = true;
    % Work to do
    irf.log('warning', ['Old intervals total        : ' num2str(nOldIntervals)])
    irf.log('warning', ['Old intervals to remove    : ' num2str(nOldObsoleteIntervals)])
    irf.log('warning', ['Old intervals to update    : ' num2str(nOldToUpdateInterval)])
    irf.log('warning', ['New intervals              : ' num2str(sum(indNewIntervals)-nOldToUpdateInterval)])
    irf.log('warning', ['Total intervals to download: ' num2str(sum(indNewIntervals))])
  end
else
  irf.log('warning', ['Dataset ' dataSet ' directory does not exist.'])
  irf.log('warning', ['Total intervals to download: ' num2str(numel(TTRequest))])
end
for j=find(indNewIntervals)'
  TTRequest.UserData(j).Status = [];
end
for j=find(~indNewIntervals)'
  TTRequest.UserData(j).Status = 1;
end
if ~exist('indexList','var') % if indexList not give, all intervals should be downloaded
  indexList = find(indNewIntervals);
end

%% Assign work
assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
if doSimulateDownload
  return;
end
%% loop through request time table
iRequest=max(indexStart,find_first_non_processed_time_interval(TTRequest));
if ~exist('indexList','var')
  nRequest	= numel(TTRequest)-iRequest+1;
  indexList	= iRequest:numel(TTRequest);
else
  nRequest	= numel(indexList);
end
while 1
  while 1 % submit next unsubmitted job
    if isempty(indexList), break; end % no more jobs
    if n_submitted_jobs(TTRequest)>=maxSubmittedJobs
      irf.log('warning','Maximum allowed number of submitted jobs is reached.');
      irf.log('warning','Pausing 10s');
      pause(10);
      break; % max allowed submitted jobs reached
    else
      % define the next job to submit
      iRequest = indexList(1);
      indexList(1) = [];
    end
    if ~isfield(TTRequest.UserData(iRequest),'Status') || ...
        ~isfield(TTRequest.UserData(iRequest),'Downloadfile') || ...
        isempty(TTRequest.UserData(iRequest).Status) || ...
        TTRequest.UserData(iRequest).Status==-1 % request not yet submitted or processed or did not succeed before
      tint=TTRequest.TimeInterval(iRequest,:);
      tint(2) = tint(2) - 1e-5; % the end is an upper boundary, to avoid the data point being in two intervals as the start and the end point we remove 10^-5 from the end time
      dataSet = TTRequest.UserData(iRequest).dataset;

      irf.log('warning',['Requesting ' dataSet ' interval #' num2str(iRequest) ...
        '(' num2str(nRequest-numel(indexList)) '/' num2str(nRequest) '): ' ...
        irf_time(tint,'tint>utc')]);
      try
        if streamData
          [download_status,downloadfile]=caa_download(tint,dataSet,...
            'stream',['downloadDirectory=' dataDir],inputParamCaaDownload{:});
        else
          [download_status,downloadfile]=caa_download(tint,dataSet,...
            'schedule','nolog','nowildcard',['downloadDirectory=' dataDir],inputParamCaaDownload{:});
        end
      catch
        download_status = -1; % something wrong with internet
        irf.log('notice','**** caa_download() DID NOT SUCCEED! ****');
      end

      if download_status == 0 % scheduling succeeded
        TTRequest.UserData(iRequest).Status=0;
        TTRequest.UserData(iRequest).Downloadfile=downloadfile;
        TTRequest.UserData(iRequest).TimeOfRequest=now;
        TTRequest.UserData(iRequest).TimeOfDownload=now;
        TTRequest.UserData(iRequest).NumberOfAttemptsToDownload=0;
        break
      elseif download_status == -1
        TTRequest.UserData(iRequest).Status=-1;
        TTRequest.UserData(iRequest).Downloadfile=[];
        TTRequest.UserData(iRequest).TimeOfRequest=now;
      end
    end
  end
  while 1 % check submitted jobs
    irf.log('notice',['Checking downloads. ' num2str(n_submitted_jobs(TTRequest)) ' jobs submitted.']);
    if n_submitted_jobs(TTRequest) == 0
      irf.log('notice','No more submitted jobs');
      break;
    end
    iSubmitted=find_first_submitted_time_interval(TTRequest);
    if isempty(iSubmitted)
      irf.log('notice','No more submitted jobs');
      break;
    end
    timeSinceDownloadSec	= (now-TTRequest.UserData(iSubmitted).TimeOfDownload)*24*3600;
    numberOfAttempts		= TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload;
    waitTimeSec				= floor(numberOfAttempts^1.3*10 - timeSinceDownloadSec) ; % wait time before downloading the file as function of number of requests already
    irf.log('warning',['Interval #' num2str(iSubmitted) ', attempt #' num2str(TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload+1)]);
    if waitTimeSec>0 ...
        && (isempty(indexList) || n_submitted_jobs(TTRequest) >= maxSubmittedJobs)
      irf.log('notice',['pause ' num2str(waitTimeSec) ' s']);
      pause(waitTimeSec);
    end
    if waitTimeSec > 0 % Avoid requesting file too fast
      break;
    end
    try
      irf.log('warning',['Interval #' num2str(iSubmitted) ': ' TTRequest.UserData(iSubmitted).Downloadfile]);
      download_status=caa_download(TTRequest.UserData(iSubmitted).Downloadfile,'nolog',['downloadDirectory=' dataDir]);
      TTRequest.UserData(iSubmitted).TimeOfDownload=now;
    catch
      download_status = -1; % something wrong with internet
      irf.log('warning','**** DID NOT SUCCEED! ****');
    end
    if download_status==1
      TTRequest.UserData(iSubmitted).Status=1; % submitted > downloaded
      irf.log('warning',['Jobs downloaded so far: ' num2str(n_downloaded_jobs(TTRequest))]);
      if mod(n_downloaded_jobs(TTRequest),10)==0 || ...%save after every 10th request
          n_submitted_jobs(TTRequest)==0 % save after last job
        varName=['TT_' dataset_mat_name(TTRequest.UserData(iSubmitted).dataset) ];
        irf.log('notice',['Saving ' varName ' to matCaaRequests/' varName]);
        eval([varName '= TTRequest;']);
        dirName='matCaaRequests';
        if ~exist(dirName,'dir'), mkdir(dirName);end
        save([dirName '/' varName],'-v7',varName);
      end
    else
      TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload=TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload+1;
      if TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload > maxNumberOfAttempts
        irf.log('warning',['*** Request #' num2str(iSubmitted) ' failed, more than 10 attempts to download ***']);
        TTRequest.UserData(iSubmitted).Status=-1;
        TTRequest.UserData(iSubmitted).Downloadfile=[];
        TTRequest.UserData(iSubmitted).TimeOfRequest=now;
        continue;
      end
      if n_submitted_jobs(TTRequest) > maxSubmittedJobs
        irf.log('warning',['More than ' num2str(maxSubmittedJobs) ' jobs sumitted, none ready, waiting 1min ....']);
        pause(60)
      end
    end
    assignin('base','TTRequest',TTRequest);
    if n_submitted_jobs(TTRequest) < maxSubmittedJobs % check next request only if less than max allowed submitted
      break
    end
  end % checking submitted jobs end
  if n_submitted_jobs(TTRequest)==0 && isempty(indexList)  %no more jobs in queue
    break;
  end
end % going through all requests
%% assign output
if nargout==1, out=TTRequest;end
%% index dataSet
local.c_update(dataSet,'datadirectory',dataDir);
%% send email when finnsihed
if sendEmailWhenFinished
  sendEmailTxt = ['local.caa_download: getting ' dataSet ' is ready ;)'];
  setpref('Internet','E_mail',sendEmailTo);
  setpref('Internet','SMTP_Server',sendEmailSmtp);
  sendmail(sendEmailTo,sendEmailTxt);
end
end
function remove_datafiles(TT,iIndex,dataDir)
% remove data files of dataset in request TT and indices iIndex (logical
% matrix)
if isa(TT,'irf.TimeTable') && islogical(iIndex)
  if numel(TT)==0
    irf.log('warning','No time intervals in request');
    return;
  elseif ~any(iIndex)
    return;
  end
  % marking files for deletion
  for j=find(iIndex)'
    fileToDelete=[dataDir filesep TT.UserData(j).filename];
    irf.log('warning',['Deleting #' num2str(j) ': ' fileToDelete]);
    system(['mv ' fileToDelete ' ' fileToDelete '.delme']);
  end
end
end
function i=find_first_non_processed_time_interval(TT)
ud=TT.UserData;
i=1;
if isfield(ud,'Status')
  while i < numel(TT) + 1
    if isempty(ud(i).Status), break; end
    if ud(i).Status == 0
      if ~isfield(ud,'Downloadfile')
        break;
      elseif isempty(ud(i).Downloadfile)
        break;
      end
    end
    i = i + 1;
  end
end
end
function i=find_first_submitted_time_interval(TT)
ud=TT.UserData;
i=[]; % default empty return
if isfield(ud,'Status')
  for j=1:numel(ud)
    if ud(j).Status==0
      i=j;
      return;
    end
  end
else
  i=1;
end
end
function n=n_submitted_jobs(TT)
if isfield(TT.UserData,'Status')
  n = sum([TT.UserData(:).Status]==0);
else
  n=0;
end
end
function n=n_downloaded_jobs(TT)
if isfield(TT.UserData,'Status')
  n = sum([TT.UserData(:).Status]==1);
else
  n=0;
end
end
function datasetMatName = dataset_mat_name(dataset)
datasetMatName = strrep(dataset,'CIS-','CIS_');
end
