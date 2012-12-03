function out=c_caa_download(varargin)
% LOCAL.C_CAA_DOWNLOAD download files to /data/caa
%
%   LOCAL.C_CAA_DOWNLOAD(dataset) download all dataset
%
% 	See also CAA_DOWNLOAD, IRF_FUNCTION_B.
%

% $Id$

if exist('/data/caa','dir')
	disp('!!!! Changing directory to /data/caa !!!');
	cd('/data/caa');
end

if nargin==1 && ischar(varargin{1})
	dataset=varargin{1};
	irf_log('dsrc','Checking list of available times');
	tt=caa_download(['list:' dataset]);
	if numel(tt)==0,
		disp('Dataset does not exist or there are no data');
		return;
	else
		irf_log('dsrc',['Checking inventory: ' irf_time(tt.TimeInterval(1,:),'tint2iso')]);
		ttInventory = caa_download(tt.TimeInterval(1,:),['list:' dataset]);
	end
	
	iData=find([ttInventory.UserData(:).number]);
	TT=select(ttInventory,iData);
	TTRequest=TT;
elseif nargin == 1 && isa(varargin{1},'TimeTable')
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
while iRequest <= numel(TTRequest),	
	tint=TTRequest.TimeInterval(iRequest,:);
	irf_log('fcal',['Requesting interval ' num2str(iRequest) '/' num2str(nRequest) ': ' irf_time(tint,'tint2iso')]);
	[download_status,downloadfile]=caa_download(tint,TTRequest.UserData(iRequest).dataset,'schedule','nolog');
	if download_status == 0, % scheduling succeeded
		TTRequest.UserData(iRequest).Status=0;
		TTRequest.UserData(iRequest).Downloadfile=downloadfile;
		TTRequest.UserData(iRequest).TimeOfRequest=now;
	end
	iRequest=iRequest+1;
	irf_log('dsrc',['=== Number of submitted jobs: ' num2str(n_submitted_jobs(TTRequest))]);
	while n_submitted_jobs(TTRequest) > 10 % number of submitted larger than 30
		irf_log('dsrc','Checking downloads...');
		iSubmitted=find_first_submitted_time_interval(TTRequest);
		irf_log('dsrc',['File: ' TTRequest.UserData(iSubmitted).Downloadfile]);
		download_status=caa_download(TTRequest.UserData(iSubmitted).Downloadfile,'nolog');
		if download_status==1,
			TTRequest.UserData(iSubmitted).Status=1; % submitted > downloaded 
			TTRequest.UserData(iSubmitted).TimeOfDownload=now;
			% save TTRequest each 10th successfull download
			irf_log('dsrc',['Jobs downloaded so far: ' num2str(n_downloaded_jobs(TTRequest))]);
			if mod(n_downloaded_jobs(TTRequest),10)==0
				irf_log('drsc',['Saving TT_' dataset ' to matCaaRequests']);
				varName=['TT_' dataset ];
				eval([varName '= TTRequest;']);
				save('CAA/matCaaRequests',varName,'-append');
			end
		else
			irf_log('proc','Waiting 1min ....');
			pause(60)
		end
	end
end

out=TTRequest;

function i=find_first_non_processed_time_interval(TT)
ud=TT.UserData;
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

