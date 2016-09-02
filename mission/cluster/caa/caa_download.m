function [downloadStatus,downloadFile]=caa_download(tint,dataset,varargin)
% CAA_DOWNLOAD Download CSA data in CDF format
%       CAA_DOWNLOAD - check the status of jobs in current directory
%
%       CAA_DOWNLOAD('list')    - list all datasets and their available times
%       TT=CAA_DOWNLOAD('list')- return time table with datasets and available times
%       CAA_DOWNLOAD('listdesc')- same with dataset description
%       CAA_DOWNLOAD('list:dataset')- list:/listdesc:  filter datasets 'dataset'
%       TT=CAA_DOWNLOAD('inventory:dataset') - return timetable with intervals when dataset has data
%
%       CAA_DOWNLOAD(tint,dataset) - download datasets matching 'dataset'
%       CAA_DOWNLOAD(tint,dataset,flags) - see different flags below
%
%       TT=CAA_DOWNLOAD(tint,'inventory') - inventory of all datasets available, return time table TT
%       TT=CAA_DOWNLOAD(tint,'inventory:dataset') - inventory of datasets matching 'dataset'
%       TT=CAA_DOWNLOAD(tint,'fileinventory:dataset') - file inventory of datasets matching 'dataset'
%
%       downloadStatus=CAA_DOWNLOAD(tint,dataset) - returns 1 if sucessfull download
%				returns 0 if request is put in the queue,
%				the information of queued requests is saved in file ".caa"
%		[downloadStatus,downloadFile]=CAA_DOWNLOAD(tint,dataset) - returns also
%				zip file link if request put in queue (good for batch processing)
%       downloadStatus=CAA_DOWNLOAD(tint,dataset,input_flags) see list of Input flags
%
%       CAA_DOWNLOAD(url_string) - download CAA data zip file from the link "url_string"
%
% Downloads CAA data in CDF format into subdirectory "CAA/"
%
%   tint   - time interval in epoch  [tint_start tint_stop]
%            or in UTC format, ex. '2005-01-01T05:00:00.000Z/2005-01-01T05:10:00.000Z'
%  dataset - dataset name, can uses also wildcard * (? is changed to *)
%
% Input flags
%   'file_interval=..' - see command line manual http://goo.gl/VkkoI, default 'file_interval=72hours'
%   'format=..'		- see command line manual http://goo.gl/VkkoI, default 'format=cdf'
%   'nowildcard'	- download the dataset without any expansion in the name and not checking if data are there
%   'overwrite'		- overwrite files in directory (to keep single cdf file)
%   'schedule'		- schedule the download, (returns zip file link)
%	'notify'        - notify by email when scheduled work is ready
%						check the readiness by executing CAA_DOWNLOAD from the same direcotry
%   'nolog'			- [default] do not log into .caa file (good for batch processing)
%   'log'			- do log into .caa file (more for interactive work)
%   'downloadDirectory=..'	- define directory for downloaded datasets (instead of default 'CAA/')
%   '&USERNAME=csaUser&PASSWORD=csaPassword'	- load data from csa using username 'uuu' and password 'ppp'
%	'cdf'           - alias of 'format=cdf'
%	'cef'           - alias of 'format=cef'
%	'json'			- return csa query in JSON format
%	'csv'			- return csa query in CSV format
%	'votable'		- return csa query in VOTABLE format
%   'stream'        - donwload using streaming interface (gzipped cef file)
%	'ingestedsince=YYYYMMDD' - download only data ingested since YYYY MM DD
%	'test'		    - test downloading some example datasets
%
%  To store your CSA user & password as defaults (e.g. 'uuu'/'ppp'):
%		datastore('csa','user','uuu')
%		datastore('csa','pwd','ppp')
%
%  Examples:
%   caa_download(tint,'list:*')       % list everything available from all sc
%   caa_download(tint,'list:*FGM*')
%   caa_download('2005-01-01T05:00:00.000Z/2005-01-01T05:10:00.000Z','list:*FGM*')
%   caa_download(tint,'C3_CP_FGM_5VPS')
%   caa_download(tint,'C?_CP_FGM_5VPS')   %    download all satellites
%
% The example list of datasets: (see also http://bit.ly/pKWVKh)
% FGM
%   caa_download(tint,'C?_CP_FGM_5VPS');
%   caa_download(tint,'C?_CP_FGM_FULL');
% EFW (L2 - full resolution, L3 - spin resolution)
%   caa_download(tint,'C?_CP_EFW_L?_E3D_INERT'); % Ex,Ey,Ez in ISR2
%   caa_download(tint,'C?_CP_EFW_L?_E3D_GSE'); % Ex,Ey,Ez in GSE
%   caa_download(tint,'C?_CP_EFW_L?_E'); % Ex,Ey in ISR2
%   caa_download(tint,'C?_CP_EFW_L?_P'); % satellite potential
%   caa_download(tint,'C?_CP_EFW_L?_V3D_GSE'); % ExB velocity GSE
%   caa_download(tint,'C?_CP_EFW_L?_V3D_INERT'); % ExB velocity ISR2
% STAFF
%   caa_download(tint,'C?_CP_STA_PSD');
%   caa_download(tint,'*STA_SM*');           % STAFF spectral matrix
% WHISPER
%   caa_download(tint,'C?_CP_WHI_NATURAL');
% CIS
%   caa_download(tint,'C?_CP_CIS_HIA_ONBOARD_MOMENTS');
%   caa_download(tint,'C?_CP_CIS_CODIF_HS_H1_MOMENTS');
%   caa_download(tint,'C?_CP_CIS_HIA_HS_1D_PEF');
%   caa_download(tint,'C?_CP_CIS_CODIF_H1_1D_PEF');
% PEACE
%   caa_download(tint,'C?_CP_PEA_PITCH_SPIN_DPFlux'); % DPFlux/DEFLux/PSD
%   caa_download(tint,'C?_CP_PEA_3DR?_PSD');
%   caa_download(tint,'C?_CP_PEA_MOMENTS')
% RAPID
%   caa_download(tint,'C?_CP_RAP_ESPCT6'); % electron omni-directional
%   caa_download(tint,'C?_CP_RAP_L3DD');   % electron, 3D distribution (standard)
%   caa_download(tint,'C?_CP_RAP_E3DD');   % electron, 3D distr. (best) in BM
%   caa_download(tint,'C?_CP_RAP_HSPCT');  % ion, omni-directional
% EPHEMERIS
%   caa_download(tint,'C?_CP_AUX_POSGSE_1M');  % position & velocity for each sc
%   caa_download(tint,'CL_SP_AUX');            % position,attitude.. for all sc
%   caa_download(tint,'C?_CP_AUX_SPIN_TIME');  % spin period, sun pulse time,..
%   caa_download(tint,'C?_JP_PMP');            % invariant latitude, MLT, L shell.


%% Check if latest irfu-matlab
% The check is appropriate to make when scientist is downloading data from CAA
persistent usingLatestIrfuMatlab

if isempty(usingLatestIrfuMatlab), % check only once if using NASA cdf
	usingLatestIrfuMatlab=irf('check');
end

%% Server defaults

% CSA
% CSA Archive Inter-Operability (AIO) System User's Manual:
% https://csa.esac.esa.int/csa/aio/html/CsaAIOUsersManual.pdf
Default.Csa.urlServer           = 'https://csa.esac.esa.int/csa/aio/';
Default.Csa.urlQuery            = 'product-action?&NON_BROWSER';
Default.Csa.urlQueryAsync       = 'async-product-action?&NON_BROWSER';
Default.Csa.urlStream           = 'streaming-action?&NON_BROWSER&gzip=1';
Default.Csa.urlInventory        = 'metadata-action?&NON_BROWSER&SELECTED_FIELDS=DATASET_INVENTORY&RESOURCE_CLASS=DATASET_INVENTORY';
Default.Csa.urlFileInventory    = 'metadata-action?&NON_BROWSER&SELECTED_FIELDS=FILE.LOGICAL_FILE_ID,FILE.START_DATE,FILE.END_DATE,FILE.CAA_INGESTION_DATE&FILE.ACTIVE=1&RESOURCE_CLASS=FILE';
Default.Csa.urlListDataset      = 'metadata-action?&NON_BROWSER&SELECTED_FIELDS=DATASET.DATASET_ID,DATASET.START_DATE,DATASET.END_DATE,DATASET.TITLE&RESOURCE_CLASS=DATASET';
Default.Csa.urlListDatasetDesc  = 'metadata-action?&NON_BROWSER&SELECTED_FIELDS=DATASET.DATASET_ID,DATASET.START_DATE,DATASET.END_DATE,DATASET.TITLE,DATASET.DESCRIPTION&RESOURCE_CLASS=DATASET';
Default.Csa.urlListFormat       = '&RETURN_TYPE=CSV';
Default.Csa.urlNotifyOn         = '';
Default.Csa.urlNotifyOff        = '&NO_NOTIFY';
Default.Csa.urlDataset          = '&DATASET_ID=';
Default.Csa.urlDataFormat       = '&DELIVERY_FORMAT=CDF';
Default.Csa.urlFileInterval     = '&DELIVERY_INTERVAL=ALL';

%% Defaults that can be overwritten by input parameters
checkDownloadStatus     = false;
checkDataInventory      = true;   % check if there are any data at caa
doLog                   = false;  % do not log into .caa file
doDataStreaming         = false;  % data streaming is in beta and supports only one dataset
doDownloadScheduling    = false;  % default download directly files
doNotifyByEmail         = false;  % default, do not notify
expandWildcards         = true;   % default is to use wildcard
overwritePreviousData   = false;  % continue adding cdf files to CAA directory
specifiedTimeInterval   = false;
specifiedFileLink       = false;
specifiedIngestedSince  = false;
downloadDirectory       = './CAA/';% local directory where to put downloaded data

%% check input
if nargin==0, % [..]=caa_download
	checkDownloadStatus=true; 
elseif  nargin == 1 && ischar(tint) && strcmpi('test',tint),  % [..]=caa_download('testcsa')
	downloadStatus = test_csa;
	if nargout == 0, clear downloadStatus;end
	return;
elseif nargout>0 && nargin>0, % [..]=caa_download(..)
	checkDownloadStatus=false;
	downloadStatus = []; % default
	doLog = false;
end

if nargin>=1, % check if first argument is not caa zip or tar.gz file link
	if ischar(tint) && ...
			(any(regexp(tint,'\.zip')) || any(regexp(tint,'\.tar.gz'))) % tint is file link
		specifiedFileLink = true;
		if nargin > 1
			varargin=[dataset varargin];
		end
	elseif ischar(tint) && any(irf_time(tint,'utc>tint')) % tint is UTC string
	elseif ischar(tint) % tint is dataset
		if nargin > 1
			varargin=[dataset varargin];
		end
		dataset=tint;
		tint=[];
	elseif isa(tint,'GenericTimeArray') && length(tint)>=2
		tintEpochUnix = [tint.start.epochUnix tint.stop.epochUnix];
		tint = tintEpochUnix;
	elseif ~isnumeric(tint)
		errStr = 'tint format not recognized';
		irf.log('critical',errStr);
		error('caa_download:tint:tint_not_defined',errStr);
	end
end
if ~isempty(varargin), % check for additional flags
	for iFlag=1:numel(varargin)
		flag=varargin{iFlag};
		if strcmpi(flag,'nowildcard'),
			expandWildcards    = false;
			checkDataInventory = false;
			doNotifyByEmail    = false;
		elseif strcmpi(flag,'overwrite'),
			overwritePreviousData = true;
		elseif any(strfind(flag,'file_interval'))
			urlFileInterval = url_parameter(flag);
		elseif any(strfind(flag,'DELIVERY_INTERVAL'))
			urlFileInterval = url_parameter(flag);
		elseif any(strfind(flag,'format'))
			urlDataFormat = url_parameter(flag);
		elseif strcmpi('cdf',flag)
			urlDataFormat = url_parameter('format=cdf');
		elseif strcmpi('cef',flag)
			urlDataFormat = url_parameter('format=cef');
		elseif any(strcmpi('schedule',flag))
			doDownloadScheduling = true;
		elseif (any(strfind(flag,'uname=')) || any(strfind(flag,'USERNAME=')) )
			urlIdentity = flag;
		elseif strcmpi('nolog',flag)
			doLog = false;
		elseif strcmpi('log',flag)
			doLog = true;
		elseif strcmpi('notify',flag)
			doNotifyByEmail = true;
		elseif any(strcmpi('csa',flag)) % download from CSA instead of CAA
			irf.log('warning','caa_download(): flag ''csa'' is not needed anymore, only CSA download available');
		elseif any(strcmpi('caa',flag)) % download from CAA instead of CAA
			irf.log('warning','caa_download(): flag ''caa'' is not supported anymore. CAA interface is closed and only CSA download is available');
		elseif strcmpi('json',flag) % set query format to JSON
			urlListFormat = '&RETURN_TYPE=JSON';
		elseif strcmpi('csv',flag) % set query format to CSV
			urlListFormat = '&RETURN_TYPE=CSV';
		elseif strcmpi('votable',flag) % set query format to VOTABLE
			urlListFormat = '&RETURN_TYPE=votable';
		elseif strfind(lower(flag),'downloaddirectory=')
			downloadDirectory = flag(strfind(flag,'=')+1:end);
			if downloadDirectory(end) ~= filesep...
					|| ~strcmp(downloadDirectory(end),'/'),
				downloadDirectory(end+1) = filesep; %#ok<AGROW>
			end
		elseif strfind(lower(flag),'ingestedsince=')
			specifiedIngestedSince = true;
			ingestedSinceYYYYMMDD = flag(strfind(flag,'=')+1:end);
		elseif any(strcmpi('stream',flag)) % data streaming
			checkDataInventory = false;
			expandWildcards		= false;
			doDataStreaming		= true;
		else
			irf.log('critical',['Flag ''' flag ''' not recognized']);
		end
	end
end

%% if needed load .caa file with status for all downloads
if doLog,
	if ~exist('.caa','file'),
		caa=cell(0);
		save -mat .caa caa;
	end
	load -mat .caa caa
end

% caa.url - links to download
% caa.dataset - dataset to download
% caa.tintiso - time interval
% caa.zip - zip files to download
% caa.status - status ('submitted','downloaded','finnished')
% caa.timeofrequest - in matlab time units

%%
if specifiedFileLink,
	if doLog % add to the submission list in caa log
		j=numel(caa)+1;
		caa{j}.url		= '*';
		caa{j}.dataset	= '*';
		caa{j}.tintiso	= '*';
		caa{j}.zip		= tint;
		caa{j}.status	= 'submitted';
		caa{j}.timeofrequest= now;
		checkDownloadStatus = true;
	else % download directly the file
		zipFileLink=tint;
		isJobFinished=get_zip_file(zipFileLink);
		if ~isJobFinished, %
			irf.log('warning','Job still not finished');
		end
		downloadStatus = isJobFinished;
		if nargout == 0,
			if downloadStatus
				disp('success!');
			else
				disp('failed!');
			end
			clear downloadStatus;
		end
		return;
	end
end

Caa = Default.Csa;
if ~exist('urlIdentity','var') || isempty(urlIdentity) % if not set by input parameters use default
	urlIdentity = get_url_identity;
end
% if variables not changed by input parameters, set them to default values
for varName = {'urlListFormat','urlDataFormat','urlFileInterval'}
	if ~exist(varName{1},'var') % if not set by input parameters use default
		eval([varName{1} ' = Caa.' varName{1} ';']);
	end
end
if doNotifyByEmail
	urlNonotify = Caa.urlNotifyOn;
else
	urlNonotify = Caa.urlNotifyOff;
end
if doDownloadScheduling
	urlQuery = Caa.urlQueryAsync;
else
	urlQuery = Caa.urlQuery;
end
if specifiedIngestedSince
	YYYY  = str2double(ingestedSinceYYYYMMDD(1:4));
	MM    = str2double(ingestedSinceYYYYMMDD(5:6));
	DD    = str2double(ingestedSinceYYYYMMDD(7:8));
	tIngestedSince = irf_time([YYYY MM DD 0 0 0]);
	urlIngestedSince = ...
		['&INGESTED_SINCE=' irf_time(tIngestedSince,'epoch>utc_yyyy-mm-ddTHH:MM:ssZ')];
else
	urlIngestedSince = '';
end
if any(strfind(urlDataFormat,'&format')),% change/add defaults, hasn't added these to above flag checking
	urlDataFormat = ['&DELIVERY_' upper(urlDataFormat(2:end))];
end
caaQuery            = [Caa.urlServer urlQuery urlIdentity urlDataFormat...
	                     urlFileInterval urlNonotify urlIngestedSince];
caaStream           = [Caa.urlServer Caa.urlStream urlIdentity urlIngestedSince];
caaInventory        = [Caa.urlServer Caa.urlInventory       urlListFormat ];
caaFileInventory    = [Caa.urlServer Caa.urlFileInventory   urlListFormat ];
caaListDataset	    = [Caa.urlServer Caa.urlListDataset     urlListFormat ];
caaListDatasetDesc  = [Caa.urlServer Caa.urlListDatasetDesc urlListFormat ];

%% Check status of downloads if needed
if checkDownloadStatus,    % check/show status of downloads from .caa file
	disp('=== status of jobs (saved in file .caa) ====');
	if ~exist('.caa','file'),
		disp('No active downloads');
		if nargout==1, downloadStatus=1; end
		return;
	else
		load -mat .caa caa
	end
	if ~isempty(caa),
		for j=1:length(caa), % go through jobs
			disp([num2str(j) '.' caa{j}.status ' ' caa{j}.dataset '-' caa{j}.tintiso]);
		end
	else
		disp('No active downloads');
		if nargout==1, downloadStatus=1; end
		return;
	end
	jobsToRemove = false(1,length(caa));
	jobsFinished = false(1,length(caa));
	for j=1:length(caa), % go through jobs
		if strcmpi(caa{j}.status,'downloaded') || strcmpi(caa{j}.status,'finnished') || strcmpi(caa{j}.status,'finished') % 'finnished shoudl be removed after some time % do nothing
			jobsFinished(j) = true;
		elseif strcmpi(caa{j}.status,'submitted'),
			disp(['=== Checking status of job nr: ' num2str(j) '==='])
			isJobFinished=get_zip_file(caa{j}.zip);
			if isJobFinished, %
				caa{j}.status='FINISHED';
				save -mat .caa caa; % changes in caa saved
			else % job not finished
				disp(['STILL WAITING TO FINISH, submitted ' num2str((now-caa{j}.timeofrequest)*24*60,3) 'min ago.']);
				if now-caa{j}.timeofrequest>1, % waiting more than 1 day
					y=input('Waiting more than 24h. Delete from list? y/n :','s');
					if strcmpi(y,'y'),
						jobsToRemove(j)=1;
					end
				end
			end
		else
			disp('ERROR: Unknown status!')
			return
		end
	end
	if sum(jobsFinished)>5, % ask for cleanup
		y=input('Shall I remove FINISHED from the list? y/n :','s');
		if strcmpi(y,'y'),
			jobsToRemove = jobsToRemove | jobsFinished;
		end
	end
	caa(jobsToRemove)=[];
	save -mat .caa caa;
	return;
end

%% check if time interval specified and define queryTime and queryTimeInventory
if isnumeric(tint) && (size(tint,2)==2), % assume tint is 2 column epoch
	tintUTC=irf_time(tint,'tint>utc_yyyy-mm-ddTHH:MM:SSZ');
	specifiedTimeInterval = true;
elseif isa(tint,'GenericTimeArray') && length(tint)==2
	tintUTC=irf_time(tint,'tint>utc_yyyy-mm-ddTHH:MM:SSZ');
	specifiedTimeInterval = true;
elseif ischar(tint), % tint is in UTC format
	tintUTC=tint;
	specifiedTimeInterval = true;
elseif isempty(tint) % will only list products
else
	disp(tint);
	error('caa_download: unknown tint format');
end

if specifiedTimeInterval
	divider=strfind(tintUTC,'/');
	t1UTC = tintUTC(1:divider-1);
	t2UTC = tintUTC(divider+1:end);
	queryTime = ['&START_DATE=' t1UTC '&END_DATE=' t2UTC];
	queryTimeFileInventory = [' AND FILE.START_DATE <= ''' t2UTC '''',...
		' AND FILE.END_DATE >= ''' t1UTC ''''];
	queryTimeInventory = [' AND DATASET_INVENTORY.START_TIME <= ''' t2UTC '''',...
		' AND DATASET_INVENTORY.END_TIME >= ''' t1UTC ''''];
end

%% define queryDataset and queryDatasetInventory
[queryDataset,queryDatasetInventory,filter] = query_dataset;

%% list data if required
if any(strfind(dataset,'list')) || any(strfind(dataset,'inventory')),     % list files
	if any(strfind(dataset,'inventory')) && ~specifiedTimeInterval
		ttTemp = caa_download(['list:' filter]);
		if isempty(ttTemp) % no dataset found
			errStr = ['Dataset ' filter ' does not exist!'];
			irf.log('critical',errStr);
			error('caa_download:dataset:doesnotexist',errStr);
		end
		if isempty(ttTemp.TimeInterval), % no time intervals to download
			irf.log('warning','No datasets to download');
			downloadStatus = ttTemp;
			return;
		end
		tint = [min(ttTemp.TimeInterval(:)) max(ttTemp.TimeInterval(:))];
		if any(strfind(dataset,'fileinventory'))
			ttTemp = caa_download(tint,['fileinventory:' filter]);
		else
			ttTemp = caa_download(tint,['inventory:' filter]);
		end
		if isfield(ttTemp.UserData(1),'number')
			iData=find([ttTemp.UserData(:).number]);
			downloadStatus=select(ttTemp,iData);
		else
			downloadStatus=ttTemp;
		end
		return
	end
	if any(strfind(dataset,'listdesc'))	% get also description
		urlListDatasets=[caaListDatasetDesc queryDatasetInventory];
		returnTimeTable='listdesc';
	elseif any(strfind(dataset,'list'))
		urlListDatasets = [caaListDataset queryDatasetInventory];
		returnTimeTable='list';
	elseif any(strfind(dataset,'fileinventory'))
		urlListDatasets = [caaFileInventory queryDatasetInventory ];
		returnTimeTable='fileinventory';
	elseif any(strfind(dataset,'inventory'))
		urlListDatasets = [caaInventory queryDatasetInventory ];
		returnTimeTable='inventory';
	end
	if specifiedTimeInterval
		if any(strfind(dataset,'fileinventory'))
			urlListDatasets = [urlListDatasets queryTimeFileInventory];
		else
			urlListDatasets = [urlListDatasets queryTimeInventory];
		end
	end
	urlListDatasets = csa_parse_url(urlListDatasets);
	irf.log('warning',['Patience! Requesting "' dataset '" ' urlListDatasets]);
	caalog=urlread(urlListDatasets);
	if isempty(caalog), % return empty output
		downloadStatus = [];
		return
	end
	if nargout == 1,
		if isempty(returnTimeTable),
			out = textscan(caalog, '%s', 'delimiter', '\n'); % cell array with lines
			downloadStatus = out{1};
		else
			downloadStatus = construct_time_table(caalog,returnTimeTable);
		end
	else
		disp(caalog);
	end
	return;
end

%% download data
% create CAA directory if needed
if ~exist('CAA','dir'), mkdir('CAA');end

if checkDataInventory
	urlListDatasets=[ caaInventory  queryDatasetInventory queryTimeInventory];
	urlListDatasets = csa_parse_url(urlListDatasets);
	irf.log('warning','Patience! Requesting list of files.');
	irf.log('notice',['URL: ' urlListDatasets]);
	caalist=urlread(urlListDatasets);
	irf.log('debug',['returned: ' caalist]);
	if isempty(caalist) % no datasets available
		irf.log('warning','There are no data sets available!');
		return;
	end
end

if doDataStreaming
	urlLine=[caaStream queryDataset queryTime];
else
	urlLine = [caaQuery queryDataset queryTime];
	irf.log('warning','Be patient! Submitting data request ...');
end

irf.log('notice',['requesting: ' urlLine]);

[status,downloadedFile] = get_zip_file(urlLine);
if nargout>=1, downloadStatus=status;end
if nargout==2, downloadFile = ''; end % default return empty
if status == 0 && exist(downloadedFile,'file')
	irf.log('warning','gz file is not returned! ');
	fid=fopen(downloadedFile);
	while 1
		tline = fgetl(fid);
		if ~ischar(tline), break, end
		disp(tline)
		if any(strfind(tline,'https:')) && any(strfind(tline,'gz')), % CSA
			downloadFile = tline(strfind(tline,'https:'):strfind(tline,'gz')+1);
		elseif any(strfind(tline,'LOGIN_REQUESTED')), % 
			disp('**ERROR**');
			disp('Your username or password are incorrect, please double check!')
			disp('See the help of caa_download on how to updated your user credentials!');
			fclose(fid);
			return;
		end
	end
	fclose(fid);
	delete(downloadedFile);
	
	if exist('downloadFile','var'),
		irf.log('warning',['Request put in queue    : ' urlLine]);
		irf.log('warning',['When ready download from: ' downloadFile]);
		if doLog
			j=length(caa)+1;
			caa{j}.url=urlLine;
			caa{j}.dataset=dataset;
			caa{j}.tintiso=tintUTC;
			caa{j}.zip = downloadFile;
			caa{j}.status = 'SUBMITTED';
			caa{j}.timeofrequest = now;
			disp('=====');
			disp('The request has been put in queue');
			disp(['When ready data will be downloaded from: ' downloadFile]);
			disp('To check the status of jobs execute: caa_download');
			caa_log({'Request put in queue: ',urlLine,...
				'When ready download from:' downloadFile});
			save -mat .caa caa
		end
		if nargout>=1, downloadStatus=0; end	% 0 if job submitted
	else
		if doLog
			disp('!!!! Did not succeed to download !!!!!');
			caa_log('Did not succeed to download');
		end
		downloadStatus=[];
	end
end

%% Nested functions
	function [status,downloadedFile]=get_zip_file(urlLink)
		% download data file, if success status=1 and file is uncompressed and moved
		% to data directory, downloadedFile is set to empty. If there is no
		% gz- data file , status=0 and downloadedFile is set to the downloaded file.
		if  ~strfind(urlLink,'.gz');  error('urlLink is not gz file!') ; end
		
		status = 0; % default
		if doDataStreaming
			% define filename
			tempFileName   = 'delme.cef';
			datasetDirName = [downloadDirectory dataset filesep];
			if ~exist(datasetDirName,'dir'),
				irf.log('notice',['Creating directory: ' datasetDirName]);
				mkdir(datasetDirName);
			end
			tempFilePath   = [datasetDirName tempFileName];
			tempFilePathGz = [tempFilePath '.gz'];
			[urlLink, tmpGetRequest] = splitUrlLink(urlLink);
			if(isempty(tmpGetRequest))
			  [downloadedFile,isReady] = urlwrite(urlLink, tempFilePathGz);
			else
			  [downloadedFile,isReady] = urlwrite(urlLink, tempFilePathGz, ...
			    'Authentication', 'Basic', 'Get', tmpGetRequest);
			end
			if isReady,
				gunzip(tempFilePathGz);
				% find the file name
				fid   = fopen(tempFilePath); % remove .gz at the end
				tline = fgetl(fid);
				while ischar(tline)
					if strfind(tline,'FILE_NAME')
						i = strfind(tline,'"');
						fileNameCefGz = [tline(i(1)+1:i(2)-1) '.gz'];
						irf.log('debug',['CEF.gz file name: ' fileNameCefGz]);
						break;
					end
					tline = fgetl(fid);
				end
				fclose(fid);
				movefile(tempFilePathGz,[datasetDirName fileNameCefGz]);
				delete(tempFilePath); % remove gunzipped file that was used only to learn the file name, otherwise cef files are kept gzipped on disc
				
				irf.log('notice',['Downloaded: ' urlLink]);
				irf.log('notice',['into ->' datasetDirName fileNameCefGz]);
				status = 1;
			else
				irf.log('warning',['Did not succed to download: ' urlLink]);
				status = 0;
			end
			return;
		end
		
		downloadedFile = [tempname '.gz'];
		[urlLink, tmpGetRequest] = splitUrlLink(urlLink);
		if(isempty(tmpGetRequest))
		  [downloadedFile,isZipFileReady] = urlwrite(urlLink, downloadedFile);
		else
		  [downloadedFile,isZipFileReady] = urlwrite(urlLink, downloadedFile, ...
		    'Authentication', 'Basic', 'Get', tmpGetRequest);
		end
		
		if isZipFileReady, %
			irf.log('notice',['Downloaded: ' urlLink]);
			irf.log('notice',['into ->' downloadedFile]);
			caa_log({'Gz file returned for request',urlLink});
			tempDirectory=tempname;
			mkdir(tempDirectory);
			try
				gunzip(downloadedFile);
				delete(downloadedFile);
				downloadedFile = downloadedFile(1:end-3); % remove '.gz'
				%filelist=untar(downloadedFile,tempDirectory);
				if( any( strcmp(version, {'8.5.0.197613 (R2015a)','8.6.0.267246 (R2015b)','9.0.0.307022 (R2016a) Prerelease'}) ) )
					if(isunix)
						% Matlab 8.5 R2015a, Matlab 8.6 R2015b and
						% 9.0 R2016a Prerelease have bugs in Matlabs untar
						% (8.5 issues with permission, 8.6 and 9.0 with 
						% excess of 100 chars).
						% try to use system built in tar on unix machines 
						irf.log('notice','Matlab R2015a/b Release or R2016a Prerelease have bug in Matlab untar, trying built in system tar on Mac/Linux.');
						cmd = sprintf('tar -xf %s --directory %s', downloadedFile, tempDirectory);
						[status, ~] = system(cmd);
						cmd = sprintf('tar -tf %s', downloadedFile);
						[~,cmdOutput] = system(cmd);
						filelist = strsplit(cmdOutput(1:end-1)); % remove end \n
						for i = 1:numel(filelist),
							filelist{i} = [tempDirectory filesep filelist{i}];
						end
						if(status~=0)
							irf.log('critical','Failed untar using built in system command. Trying Matlab, but this may fail with long file names.');
							filelist = untar(downloadedFile, tempDirectory);
						end
					else
						irf.log('warning','Matlab R2015a/b Release or R2016a Prerelease on Windows, trying built in Matlab untar. This may fail with long file names!');
						filelist = untar(downloadedFile, tempDirectory);
					end
				else
					filelist = untar(downloadedFile, tempDirectory);
				end
				
				if isempty(filelist)
					irf.log('warning','Returned gz file is empty');
					caa_log('gz file empty.');
				else
					move_to_caa_directory(filelist);
				end
				status=1;
				delete(downloadedFile);
				downloadedFile = '';
			catch
				irf.log('critical','Invalid gz file')
			end
			rmdir(tempDirectory,'s');
		else
			irf.log('warning',['There is no gz file: ' urlLink]);
		end
	end

	% Nested function
	function move_to_caa_directory(filelist)
		for jj=1:length(filelist),
			isDataSet = ~any(strfind(filelist{jj},'log'));
			if isDataSet, % dataset files (cdf_convert_summary.log not copied)
				ii=strfind(filelist{jj},filesep);
				dataset=filelist{jj}(ii(end-1)+1:ii(end)-1);
				datasetDirName = [downloadDirectory dataset];
				if ~exist(datasetDirName,'dir'),
					irf.log('notice',['Creating directory: ' datasetDirName]);
					mkdir(datasetDirName);
				elseif overwritePreviousData
					delete([datasetDirName filesep '*']);
				end
				irf.log('notice',['file:      ' filelist{jj}]);
				irf.log('notice',['moving to directory: ' datasetDirName]);
				movefile(filelist{jj},datasetDirName);
			end
		end
	end

	% Nested function. Check url parameter syntax in string
	function paramOut=url_parameter(paramIn)
		if paramIn(1)~= '&'
			paramOut=['&' paramIn];
		else
			paramOut=paramIn;
		end;
	end

	% Nested function. Save text to log file
	function caa_log(logText)
		tt=irf_time;
		if ischar(logText), logText={logText};end
		if iscellstr(logText)
			fid=fopen('.caa.log','a');
			if fid==-1 % cannot open .caa.log
				if isempty(whos('-file','.caa','logFileName')) % no log file name in caa
					logFileName=tempname;
					save -append .caa logFileName;
				else
					load -mat .caa logFileName;
				end
				fid=fopen(logFileName,'a');
				if fid==-1,
					irf.log('critical','log file cannot be opened, no log entry');
					return;
				end
			end
			fprintf(fid,'\n[%s]\n',tt);
			for jLine=1:numel(logText),
				fprintf(fid,'%s\n',logText{jLine});
			end
			fclose(fid);
		end
	end

	% Nested function. Parse url line syntax before requesting from www
	function out = csa_parse_url(in)
		% replace all space and stars with %20 and %25 correspondingly
		out = strrep(in,'*','%25');
		out = strrep(out,' ','%20');
	end

	% Nested function.
	function [queryDataset,queryDatasetInventory,filter] = query_dataset
		% for wildcards, inventory requests use '%' as wildcard,
		% while data requests use '*' (something that was not easy to implement)
		if any(strfind(dataset,'list')) || any(strfind(dataset,'inventory')),     % list files
			if strcmpi(dataset,'list') && strcmpi(dataset,'listdesc'), % list all files
				filter='*';
			else                        % list only filtered files
				filter=dataset(strfind(dataset,':')+1:end);
			end
		else
			filter = dataset;
		end
		if expandWildcards,
			filter(strfind(filter,'?'))='*'; % substitute  ? to * (to have the same convention as in irf_ssub)
			if (any(strfind(filter,'CIS')) || any(strfind(filter,'CODIF')) || any(strfind(filter,'HIA')))
				filter(strfind(filter,'_'))='*'; % substitute  _ to * (to handle CIS products that can have - instead of _)
			end
		end
		queryDataset = ['&DATASET_ID=' filter];
		queryDatasetInventory = ['&QUERY=DATASET.DATASET_ID like ''' csa_parse_url(filter) ''''];
	end

	% Nested function. Get CSA identity
	function urlIdentity = get_url_identity
		csaUser = datastore('csa','user');
		if isempty(csaUser)
			csaUser = input('Input csa username [default:avaivads]:','s');
			if isempty(csaUser),
				disp('Please register at ______? and later use your username and password.');
				csaUser='avaivads';
			end
			datastore('csa','user',csaUser);
		end
		csaPwd = datastore('csa','pwd');
		if isempty(csaPwd)
			csaPwd = input('Input csa password [default:!kjUY88lm]:','s');
			if isempty(csaPwd), csaPwd='!kjUY88lm';end
			datastore('csa','pwd',csaPwd);
		end
		urlIdentity = ['&USERNAME=' csaUser '&PASSWORD=' csaPwd];
	end

	% Nested function. 
	function TT=construct_time_table(caalog,returnTimeTable)
		TT=irf.TimeTable;
		switch returnTimeTable
			case 'inventory'
				%"DATASET_INVENTORY.DATASET_ID","DATASET_INVENTORY.START_DATE","DATASET_INVENTORY.END_DATE"
				textLine=textscan(caalog,'"%[^"]","%[^"]","%[^"]","%[^"]","%[^"]"');
				TT.UserData(numel(textLine{1})-1).dataset = [];
				[TT.UserData(:).dataset]=deal(textLine{1}{2:end});
				for jj = 1:numel(TT.UserData),
					TT.UserData(jj).number = str2double(textLine{4}{1+jj});
				end
				[TT.UserData(:).version]=deal(textLine{5}{2:end});
			case 'fileinventory'
				%"FILE.LOGICAL_FILE_ID","FILE.START_DATE","FILE.END_DATE","FILE.CAA_INGESTION_DATE"
				textLine=textscan(caalog,'"%[^"]","%[^"]","%[^"]","%[^"]"');
				TT.UserData(numel(textLine{1})-1).dataset = [];
				[TT.UserData(:).filename]=deal(textLine{1}{2:end});
				[TT.UserData(:).caaIngestionDate]=deal(textLine{4}{2:end});
			case 'list'
				%DATASET.DATASET_ID,DATASET.START_DATE,DATASET.END_DATE,DATASET.TITLE
				textLine=textscan(caalog,'"%[^"]","%[^"]","%[^"]","%[^"]"');
				TT.UserData(numel(textLine{1})-1).dataset = [];
				[TT.UserData(:).dataset]=deal(textLine{1}{2:end});
			case 'listdesc'
				%DATASET.DATASET_ID,DATASET.START_DATE,DATASET.END_DATE,DATASET.TITLE,DATASET.DESCRIPTION
				textLine=textscan(caalog,'"%[^"]","%[^"]","%[^"]","%[^"]","%[^"]"');
				TT.UserData(numel(textLine{1})-1).dataset = [];
				[TT.UserData(:).dataset]=deal(textLine{1}{2:end});
				[TT.UserData(:).description] = deal(textLine(:).description);
			otherwise
				return;
		end
		tStart = arrayfun(@(x) irf_time(x{1},'utc>epoch'),textLine{2}(2:end));
		tEnd   = arrayfun(@(x) irf_time(x{1},'utc>epoch'),textLine{3}(2:end));
		tint = [tStart tEnd];
		TT.TimeInterval=tint;
		TT.Header = {};
		TT.Comment=cell(numel(TT),1);
		TT.Description=cell(numel(TT),1);
	end
end
%% Functions (not nested)
function ok = test_csa
ok=false;
disp('--- TESTING CSA ---');
currentDir = pwd;
c1 = onCleanup(@()cd(currentDir));
tempDir=tempname;
mkdir(tempDir);
cd(tempDir);
disp(['Moving to temporary directory: ' tempDir]);
c2 = onCleanup(@() rmdir(tempDir,'s'));

% 30s of data
tintUTC = '2005-01-01T05:00:00.000Z/2005-01-01T05:00:30.000Z';
tt = caa_download(tintUTC,'list:C3_CP_FGM*');
if isempty(tt), disp('FAILED!'); return; end
ok = caa_download(tintUTC,'C3_CP_FGM_5VPS');
if ~ok, disp('FAILED!'); return; end
ok = caa_download(tintUTC,'C3_CP_FGM_5VPS','stream');
if ~ok, disp('FAILED!'); return; end
ok = caa_download(tintUTC,'C?_CP_FGM_5VPS');
if ~ok, disp('FAILED!'); return; end
disp('--- TEST PASSED! ---');
ok=true;
end

function [url, getRequest] = splitUrlLink(urlLink)
% Help function to split CSA url requests to account for new interface.
% CSA is to replace thier interface from 2016/05/04 onward.
  if(~isempty(regexpi(urlLink,'password')))
    % It it a password protected page being requested, split it.
    tmpSplit = strsplit(urlLink,'&');
    url = strrep(tmpSplit{1}, '?',''); % strip "?"
    for ii=2:length(tmpSplit)
      tmpSplit2 = strsplit(tmpSplit{ii},'=');
       % Single element, such as "NON_BROWSER", Add "1" as second argument.
      if(size(tmpSplit2,2)==1), tmpSplit2{2}='1'; end
      if(exist('getRequest','var'))
        getRequest = [getRequest, tmpSplit2]; %#ok<AGROW>
      else
        getRequest = tmpSplit2;
      end
    end
  else
    % It is not password protected URL, return unaltered urlLink (used for
    % urlread and not urlwrite?).
    url = urlLink;
    getRequest = [];
  end
end
