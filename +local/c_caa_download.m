function out=c_caa_download(varargin)
% LOCAL.C_CAA_DOWNLOAD download files to /data/caa
%
%   LOCAL.C_CAA_DOWNLOAD(dataset) download all dataset
%
% 	See also CAA_DOWNLOAD, IRF_FUNCTION_B.
%

% $Id$

maxSubmittedJobs = 13;
maxNumberOfAttempts = 20;

if exist('/data/caa','dir')
	disp('!!!! Changing directory to /data/caa !!!');
	cd('/data/caa');
end

if nargin==1 && ischar(varargin{1})
	dataSet=varargin{1};
	irf_log('dsrc','Checking list of available times');
	tt=caa_download(['list:' dataSet]);
	if numel(tt)==0,
		disp('Dataset does not exist or there are no data');
		return;
	else
		irf_log('dsrc',['Checking inventory: ' irf_time(tt.TimeInterval(1,:),'tint2iso')]);
		ttInventory = caa_download(tt.TimeInterval(1,:),['list:' dataSet]);
	end

	iData=find([ttInventory.UserData(:).number]);
	TT=select(ttInventory,iData);
	TTRequest=TT;
	assignin('base','TTRequest',TTRequest); % TTRequest assign so that one can work
elseif nargin == 1 && isa(varargin{1},'irf.TimeTable')
	TTRequest=varargin{1};
else
	irf_log('fcal','See syntax: help local.c_caa_download');
	return;
end


% TTRequest.UserData.Status = 1 - downloaded, 0 - submitted, empty - not processed
% TTRequest.UserData.Downloadfile zip file to download (important if status=0)
% TTRequest.UserData.TimeOfRequest
% TTRequest.UserData.TimeOfDownload
iRequest=find_first_non_processed_time_interval(TTRequest);
nRequest=numel(TTRequest)-iRequest+1;
while 1
  if iRequest <= numel(TTRequest) && n_submitted_jobs(TTRequest)<maxSubmittedJobs
	tint=TTRequest.TimeInterval(iRequest,:);
	irf_log('fcal',['Requesting interval ' num2str(iRequest) '(' num2str(nRequest-(numel(TTRequest)-iRequest)) '/' num2str(nRequest) '): ' irf_time(tint,'tint2iso')]);
	dataSet = TTRequest.UserData(iRequest).dataset;
	try 
	  [download_status,downloadfile]=caa_download(tint,dataSet,'schedule','nolog','nowildcard');
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
	  elseif download_status == -1,
		TTRequest.UserData(iRequest).Status=-1;
		TTRequest.UserData(iRequest).Downloadfile=[];
		TTRequest.UserData(iRequest).TimeOfRequest=now;
		iRequest = iRequest + 1;
		continue;
	end
	iRequest=iRequest+1;
  end
	while 1   % check submitted jobs
		irf_log('dsrc',['Checking downloads. ' num2str(n_submitted_jobs(TTRequest)) ' jobs submitted.']);
		iSubmitted=find_first_submitted_time_interval(TTRequest);
		if now-TTRequest.UserData(iSubmitted).TimeOfDownload>(1+TTRequest.UserData(iSubmitted).NumberOfAttemptsToDownload)*1/24/60 % more than 1min*"number of attempts" since last request 
		  try
			irf_log('dsrc',['File #' num2str(iSubmitted) ': ' TTRequest.UserData(iSubmitted).Downloadfile]);
			download_status=caa_download(TTRequest.UserData(iSubmitted).Downloadfile,'nolog');
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
				varName=['TT_' TTRequest.UserData(iSubmitted).dataset ];
				irf_log('drsc',['Saving ' varName ' to CAA/matCaaRequests/' varName]);
				eval([varName '= TTRequest;']);
				dirName=['CAA/matCaaRequests'];
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
	  else
		pause(10)
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

if nargout==1, out=TTRequest;end

function i=find_first_non_processed_time_interval(TT)
ud=TT.UserData;
i=numel(TT)+1; % default all is processed
if isfield(ud,'Status')
	for j=1:numel(ud),
		if isempty(ud(j).Status)
			i=j;
			return;
		end
	end
else
	i=1;
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

