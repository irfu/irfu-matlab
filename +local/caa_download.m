function out=caa_download(varargin)
% LOCAL.CAA_DOWNLOAD download full datasets from CAA
% Downloads all data from CAA database, in case data
% already exists on disk, downloads only newer version files.
% Dataset location - /data/caalocal
% Dataset request timetable information - /data/caalocal/matCaaRequests
% Index location - /data/caalocal/index
% Current request time table is accessible in variable TTRequest
%
%   LOCAL.CAA_DOWNLOAD(dataset) download all dataset
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
%

% Request time table TTRequest structure
%
% TTRequest.UserData.Status = 1 - downloaded, 0 - submitted, empty - not processed
% TTRequest.UserData.Downloadfile zip file to download (important if status=0)
% TTRequest.UserData.TimeOfRequest
% TTRequest.UserData.TimeOfDownload
% TTRequest.UserData.NumberOfAttemptsToDownload
% TTRequest.UserData.dataset
% TTRequest.UserData.number - number of entries
% TTRequest.UserData.version - version of dataset

%% Defaults
dataDirectory = '/data/caalocal';
maxSubmittedJobs = 13;
maxNumberOfAttempts = 20;
isInputDatasetName = false;
sendEmailWhenFinished = false;
streamData = false;               % download cdf files asynchronously
% so far undocumented feature
% use datastore info in local to send email when finnished
if exist('sendmail','file')==2,
	if isempty(fields(datastore('local'))),
		% nothing in datastore('local'), do not send email
	else
		sendEmailWhenFinished = true;
		sendEmailFrom = datastore('local','sendEmailFrom');
		if isempty(sendEmailFrom),
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
		if isempty(sendEmailSmtp),
			disp('Please define your local SMTP server. ');
			disp('For example in IRFU: sol.irfu.se');
			disp('Execute in matlab (adjust accordingly):')
			disp(' ');
			disp('>  datastore(''local'',''sendEmailSmtp'',''sol.irfu.se'')');
			disp(' ');
			return;
		end
		sendEmailTo = datastore('local','email');
		if isempty(sendEmailTo),
			disp('Please specify your email.');
			disp('For example: name@gmail.com');
			disp('Execute in matlab (adjust accordingly):')
			disp(' ');
			disp('>  datastore(''local'',''email'',''name@gmail.com'')');
			disp(' ');
			return;
		end
	end
end
%% change to data directory
if exist(dataDirectory,'dir')
	disp(['!!!! Changing directory to ' dataDirectory ' !!!']);
	cd(dataDirectory);
else
	disp(['DOES NOT EXIST data directory: ' dataDirectory ' !!!']);
	out=[];
	return;
end
%% check input: get inventory and construct time table if dataset
if nargin == 0,
  help local.caa_download
  return
end
if ischar(varargin{1})
	dataSet=varargin{1};
	isInputDatasetName = true;
	irf_log('dsrc','Checking list of available times');
	TT=caa_download(['listdata:' dataSet]);
	if numel(TT)==0,
		disp('Dataset does not exist or there are no data');
		return;
	end
	TTRequest=TT;
	assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
elseif isa(varargin{1},'irf.TimeTable')
	TTRequest=varargin{1};
	dataSet=dataset_mat_name(TTRequest.UserData(1).dataset);
else
	irf_log('fcal','See syntax: help local.c_caa_download');
	return;
end
if nargin >= 2 && ischar(varargin{2}) && strcmpi(varargin{2},'stream')
	streamData = true; 
	irf_log('fcal','Streaming data from CAA.');
end
%% check which time intervals are already downloaded, remove obsolete ones
requestListVariableName=['TT_' dataSet ];
requestListDirectory='matCaaRequests';
requestListVariableFile=[requestListDirectory filesep requestListVariableName '.mat'];
if ~exist(requestListDirectory,'dir'),
	mkdir(requestListDirectory);
else % merge request lists
	if exist(requestListVariableFile,'file'),
		load(requestListVariableFile,requestListVariableName);
		TTRequest_old=eval(requestListVariableName);
		irf_log('fcal','Previous request list exists, merging...');
		[~,iiobsolete]=setdiff(TTRequest_old,TTRequest); % those intervals that are not anymore in inventory
		irf_log('fcal',['Obselete intervals in old requests list: ' num2str(iiobsolete)]);
		remove_datafiles(TTRequest_old,iiobsolete);
		TTRequest_old=remove(TTRequest_old,iiobsolete);
		[~,iinew]=setdiff(TTRequest,TTRequest_old); % new time intervals
		irf_log('fcal',['New time intervals, not in old requests list: ' num2str(iinew)]);
		[~,iiold,iinew]=common(TTRequest_old,TTRequest); % common time intervals
		updateFields={'Status','TimeOfRequest','TimeOfDownload','NumberOfAttemptsToDownload'};
		for ii=1:numel(iiold),
			if TTRequest_old.UserData(iiold(ii)).version == TTRequest.UserData(iinew(ii)).version && ...
					TTRequest_old.UserData(iiold(ii)).number == TTRequest.UserData(iinew(ii)).number
				for jj=1:numel(updateFields),
					if isInputDatasetName
						TTRequest.UserData(iinew(ii)).(updateFields{jj})=TTRequest_old.UserData(iiold(ii)).(updateFields{jj});
					end
				end
			else
				remove_datafiles(TTRequest_old,iiold(ii));
				irf_log('fcal',['Request download for interval ' iinew(ii)]);
			end
		end
	end
end

assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
%% loop through request time table
iRequest=find_first_non_processed_time_interval(TTRequest);
nRequest=numel(TTRequest)-iRequest+1;
while 1
	while 1 % submit next unsubmitted job
		if iRequest > numel(TTRequest), break;end % no more jobs
		if n_submitted_jobs(TTRequest)>=maxSubmittedJobs, break; end % max allowed submitted jobs reached
		if ~isfield(TTRequest.UserData(iRequest),'Status') || ...
		  		~isfield(TTRequest.UserData(iRequest),'Downloadfile') || ...
				isempty(TTRequest.UserData(iRequest).Status) || ...
				TTRequest.UserData(iRequest).Status==-1 % request not yet submitted or processed or did not succeed before
			tint=TTRequest.TimeInterval(iRequest,:);
			irf_log('fcal',['Requesting interval ' num2str(iRequest) '(' num2str(nRequest-(numel(TTRequest)-iRequest)) '/' num2str(nRequest) '): ' irf_time(tint,'tint2iso')]);
			dataSet = TTRequest.UserData(iRequest).dataset;
			try
				if streamData
				[download_status,downloadfile]=caa_download(tint,dataSet,'stream',['downloadDirectory=' dataDirectory]);
				else
				[download_status,downloadfile]=caa_download(tint,dataSet,'schedule','nolog','nowildcard');
				end
			catch
				download_status = -1; % something wrong with internet
				irf_log('dsrc','**** DID NOT SUCCEED! ****');
			end
			if download_status == 0, % scheduling succeeded
				TTRequest.UserData(iRequest).Status=0;
				TTRequest.UserData(iRequest).Downloadfile=downloadfile;
				TTRequest.UserData(iRequest).TimeOfRequest=now;
				TTRequest.UserData(iRequest).TimeOfDownload=now;
				TTRequest.UserData(iRequest).NumberOfAttemptsToDownload=0;
				iRequest=iRequest+1;
				break
			elseif download_status == -1,
				TTRequest.UserData(iRequest).Status=-1;
				TTRequest.UserData(iRequest).Downloadfile=[];
				TTRequest.UserData(iRequest).TimeOfRequest=now;
			end
		end
		iRequest=iRequest+1;
	end
	while 1 % check submitted jobs
		irf_log('dsrc',['Checking downloads. ' num2str(n_submitted_jobs(TTRequest)) ' jobs submitted.']);
		if n_submitted_jobs(TTRequest) == 0, 
		  irf_log('fcal','No more submitted jobs');
		  break;
		end
		iSubmitted=find_first_submitted_time_interval(TTRequest);
		if isempty(iSubmitted),
			irf_log('fcal','No more submitted jobs');
			break;
		end
		timeSinceDownloadSec	= (now-TTRequest.UserData(iSubmitted).TimeOfDownload)/24/3600;
		numberOfAttempts		= TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload;
		waitTimeSec				= numberOfAttempts*60 - timeSinceDownloadSec ; % 'number of attempts' * minutes is expected wait time before next check
		irf_log('dsrc',['Job #' num2str(iSubmitted) ', attempt #' num2str(TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload+1)]);
		if waitTimeSec>0,
			irf_log('dsrc',['pause ' num2str(waitTimeSec) ' s']);
			pause(waitTimeSec);
		end
		try
			irf_log('dsrc',['File #' num2str(iSubmitted) ': ' TTRequest.UserData(iSubmitted).Downloadfile]);
			download_status=caa_download(TTRequest.UserData(iSubmitted).Downloadfile,'nolog',['downloadDirectory=' dataDirectory]);
			TTRequest.UserData(iSubmitted).TimeOfDownload=now;
		catch
			download_status = -1; % something wrong with internet
			irf_log('dsrc','**** DID NOT SUCCEED! ****');
		end
		if download_status==1,
			TTRequest.UserData(iSubmitted).Status=1; % submitted > downloaded
			irf_log('dsrc',['Jobs downloaded so far: ' num2str(n_downloaded_jobs(TTRequest))]);
			if mod(n_downloaded_jobs(TTRequest),10)==0 || ...%save after every 10th request
					n_submitted_jobs(TTRequest)==0 % save after last job
				varName=['TT_' dataset_mat_name(TTRequest.UserData(iSubmitted).dataset) ];
				irf_log('drsc',['Saving ' varName ' to matCaaRequests/' varName]);
				eval([varName '= TTRequest;']);
				dirName='matCaaRequests';
				if ~exist(dirName,'dir'), mkdir(dirName);end
				save([dirName '/' varName],'-v7',varName);
			end
		else
			TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload=TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload+1;
			if TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload > maxNumberOfAttempts
				irf_log('dsrc',['*** Request #' num2str(iSubmitted) ' failed, more than 10 attempts to download ***']);
				TTRequest.UserData(iSubmitted).Status=-1;
				TTRequest.UserData(iSubmitted).Downloadfile=[];
				TTRequest.UserData(iSubmitted).TimeOfRequest=now;
				continue;
			end
			if n_submitted_jobs(TTRequest) > maxSubmittedJobs
				irf_log('proc',['More than ' num2str(maxSubmittedJobs) ' jobs sumitted, none ready, waiting 1min ....']);
				pause(60)
			end
		end
		assignin('base','TTRequest',TTRequest);
		if n_submitted_jobs(TTRequest) < maxSubmittedJobs; % check next request only if less than max allowed submitted
			break
		end
	end % checking submitted jobs end
	if n_submitted_jobs(TTRequest)==0   %no more jobs in queue
		break;
	end
end % going through all requests
%% assign output
if nargout==1, out=TTRequest;end
%% send email when finnsihed
if sendEmailWhenFinished
	sendEmailTxt = ['local.caa_download: getting ' dataSet ' is ready ;)'];
	setpref('Internet','E_mail',sendEmailTo);
	setpref('Internet','SMTP_Server',sendEmailSmtp);
	sendmail(sendEmailTo,sendEmailTxt);
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
function i=find_first_submitted_time_interval(TT)
ud=TT.UserData;
i=[]; % default empty return
if isfield(ud,'Status')
	for j=1:numel(ud),
		if ud(j).Status==0
			i=j;
			return;
		end
	end
else
	i=1;
end
function n=n_submitted_jobs(TT)
if isfield(TT.UserData,'Status')
	n = sum([TT.UserData(:).Status]==0);
else
	n=0;
end
function n=n_downloaded_jobs(TT)
if isfield(TT.UserData,'Status')
	n = sum([TT.UserData(:).Status]==1);
else
	n=0;
end
function ok=remove_datafiles(TT,ii)
% remove data files of dataset in request TT and indices ii
ok=false; % default
if isa(TT,'irf.TimeTable') && isnumeric(ii)
	if numel(TT)==0,
		irf_log('fcal','No time intervals in request');
		return;
	elseif isempty(ii)
		ok=true;
		return;
	end
	TTremove=select(TT,ii);
	dataSet=TTremove.UserData(1).dataset;
	% get dataset file index
	load('caa',['index_' dataSet]);
	index=eval(['index_' dataSet]);
	indexTT=irf.TimeTable([index.tstart(:) index.tend(:)]);
	for kk=1:numel(indexTT)
		indexTT.UserData(kk).filename=index.filename(kk,:);
	end
	% check which files to remove
	[~,iTT,iIndex]=common(TTremove,indexTT);
	for j=1:numel(iIndex)
		fileToDelete=indexTT.UserData(iIndex(j)).filename;
		irf_log('fcal',['Deleting #' num2str(iIndex(j)) ': ' fileToDelete]);
		eval(['!mv ' fileToDelete ' ' fileToDelete '.delme']);
	end
	if numel(iTT) == numel(TTremove)
		ok=true;
	end
end
function datasetMatName = dataset_mat_name(dataset)
datasetMatName = strrep(dataset,'CIS-','CIS_');

