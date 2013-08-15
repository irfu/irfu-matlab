function [download_status,downloadfile]=caa_download(tint,dataset,varargin)
% CAA_DOWNLOAD Download CAA data in CDF format
%       CAA_DOWNLOAD - check the status of jobs in current directory
%
%       CAA_DOWNLOAD('list')    - list all datasets and their available times
%       TT=CAA_DOWNLOAD('list')- return time table with datasets and available times
%       CAA_DOWNLOAD('listdesc')- same with dataset description
%       CAA_DOWNLOAD('listgui') - same presenting output in separate window
%       CAA_DOWNLOAD('list:dataset')- list:/listdesc:/listgui:  filter datasets 'dataset'
%       TT=CAA_DOWNLOAD('listdata:dataset') - return timetable with intervals when dataset has data
%
%       CAA_DOWNLOAD(tint,dataset) - download datasets matching 'dataset'
%       CAA_DOWNLOAD(tint,dataset,flags) - see different flags below
%
%       TT=CAA_DOWNLOAD(tint,'list') - inventory of all datasets available, can return time table TT
%       TT=CAA_DOWNLOAD(tint,'list:dataset') - only inventory datasets matching 'dataset'
%
%       download_status=CAA_DOWNLOAD(tint,dataset) - returns 1 if sucessfull download
%				returns 0 if request is put in the queue,
%				the information of queued requests is saved in file ".caa"
%		[download_status,downloadfile]=CAA_DOWNLOAD(tint,dataset) - returns also
%				zip file link if request put in queue (good for batch processing)
%       download_status=CAA_DOWNLOAD(tint,dataset,input_flags) see list of Input flags
%
%       CAA_DOWNLOAD(url_string) - download CAA data zip file from the link "url_string"
%
% Downloads CAA data in CDF format into subdirectory "CAA/"
%
%   tint   - time interval in epoch  [tint_start tint_stop]
%            or in ISO format, ex. '2005-01-01T05:00:00.000Z/2005-01-01T05:10:00.000Z'
%  dataset - dataset name, can uses also wildcard * (? is changed to *)
%
% Input flags
%   'file_interval=..' - see command line manual http://goo.gl/VkkoI, default 'file_interval=72hours'
%   'format=..'		- see command line manual http://goo.gl/VkkoI, default 'format=cdf'
%   'nowildcard'	- download the dataset without any expansion in the name and not checking if data are there
%   'overwrite'		- overwrite files in directory (to keep single cdf file)
%   'schedule'		- schedule the download, (returns zip file link)
%						check the readiness by executing CAA_DOWNLOAD from the same direcotry
%   'nolog'			- do not log into .caa file (good for batch processing)
%   'downloadDirectory=..'	- define directory for downloaded datasets (instead of deaful 'CAA/')
%   'uname=uuu&pwd=ppp'	- load data from caa using username 'uuu' and password 'ppp'
%   'cfa'           - download from CFA instead of CAA, basic functions are
%                     implemented. Downloaded files go into './CFA/' folder.
%                     This must be taken into account when loading data.
%
%
%  To store your caa user & password as defaults (e.g. 'uuu'/'ppp'): 
%		datastore('caa','user','uuu')
%		datastore('caa','pwd','ppp')
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

% Test flags
%   'test'						- use caa test server instead

 
%% Check if latest irfu-matlab 
% The check is appropriate to make when scientist is downloading data from CAA
persistent usingLatestIrfuMatlab urlIdentity 
global downloadFromCFA

if isempty(usingLatestIrfuMatlab), % check only once if using NASA cdf
	usingLatestIrfuMatlab=irf('check');
end

%% Defaults
checkDownloadStatus		= false;
doLog					= true;			% log into .caa file
overwritePreviousData	= false;		% continue adding cdf files to CAA directory
expandWildcards			= true;			% default is to use wildcard
checkIfDataAreAtCaa		= true;			% check if there are any data at caa
urlNonotify				= '&nonotify=1';% default is not notify by email
urlFileInterval='&file_interval=72hours'; % default time interval of returned files
urlSchedule='';							% default do not have schedule option
urlFormat='&format=cdf';				% default is CDF (3.3) format
caaServer='http://caa.estec.esa.int/';	% default server
urlIdentity=get_url_identity('caa');	% default identity, default caa
downloadDirectory = './CAA/';			% local directory where to put downloaded data, default in current directory under 'CAA' subdirectory
downloadFromCFA = 0;                    % default is to download from CAA
%% load .caa file with status for all downloads
if doLog,
	if exist('.caa','file') == 0,
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
%% check input
if nargin==0, checkDownloadStatus=true; end
if nargout>0 && nargin>0, 
	checkDownloadStatus=false; 
	doLog = false;
end
if nargin>2, % check for additional flags
	for iFlag=1:numel(varargin)
		flag=varargin{iFlag};
		if strcmpi(flag,'test'),  % use test server
			caaServer='http://caa5.estec.esa.int/caa_query/';
			urlNonotify='';           % notify also by email
		elseif strcmpi(flag,'nowildcard'),
			expandWildcards = false;
			checkIfDataAreAtCaa = false;
			urlNonotify = '&nonotify=1';
		elseif strcmpi(flag,'overwrite'),
			overwritePreviousData = true;
		elseif any(strfind(flag,'file_interval'))
			urlFileInterval = urlparameter(flag);
		elseif any(strfind(flag,'format'))
			urlFormat = urlparameter(flag);
		elseif any(strcmpi('schedule',flag))
			urlSchedule = '&schedule=1';
		elseif any(strfind(flag,'uname='))
			urlIdentity = flag;
		elseif any(strcmpi('nolog',flag))
			doLog = false;
		elseif any(strcmpi('inventory',flag))
			urlSchedule = '&inventory=1';
        elseif any(strcmpi('cfa',flag)) % download from CFA instead of CAA
			downloadFromCFA = 1;
            urlIdentity = get_url_identity('cfa');
		elseif any(strfind(flag,'USERNAME='))
			urlIdentity = flag;
        elseif any(strcmpi('json',flag)) % set query format to 
			queryFormat = '&RETURN_TYPE=JSON';
        elseif any(strcmpi('csv',flag)) % set query format to 
			queryFormat = '&RETURN_TYPE=CSV';
        elseif any(strcmpi('votable',flag)) % set query format to 
			queryFormat = '&RETURN_TYPE=votable';
		elseif any(strfind(lower(flag),'downloaddirectory='))
			downloadDirectory = flag(strfind(flag,'=')+1:end);
			if downloadDirectory(end) ~= filesep,
				downloadDirectory(end+1) = filesep;
			end
		else
			irf_log('fcal',['Flag ''' flag ''' not recognized']);
		end
	end
end
if nargin>=1, % check if first argument is not caa zip file link
	if ischar(tint) && any(regexp(tint,'\.zip'))% tint zip file link
		if nargin>1 && ischar(dataset) && strcmpi(dataset,'nolog')
			doLog=false;
		end
		if doLog
			j=numel(caa)+1;
			caa{j}.url		= '*';
			caa{j}.dataset	= '*';
			caa{j}.tintiso	= '*';
			caa{j}.zip		= tint;
			caa{j}.status	= 'submitted';
			caa{j}.timeofrequest= now;
			checkDownloadStatus = true;
		else
			zipFileLink=tint;
			isJobFinished=get_zip_file(zipFileLink);
			if ~isJobFinished, %
				irf_log('dsrc','Job still not finished');
			end
			download_status = isJobFinished;
			return;
		end
	elseif ischar(tint) && any(irf_time(tint,'iso2tint')) % tint is tintiso 
	elseif ischar(tint) % tint is dataset 
		dataset=tint;
		tint=[];
	elseif ~isnumeric(tint)
		help caa_download;
		return;
	end
end
caaQuery=[caaServer 'caa_query/?'];
caaInventory=[caaServer 'cgi-bin/inventory.cgi/?'];

%%%%%%%%% CFA stuff, 
% CFA Archive Inter-Operability (AIO) System User's Manual: 
%   http://cfadev.esac.esa.int/cfa/aio/html/CfaAIOUsersManual.pdf
if downloadFromCFA % change/add defaults, hasn't ad ded these to above flag checking
    % change
    urlNonotify = '&NO_NOTIFY';                % default is not notify by email        
    urlFormat = ['&DELIVERY_' upper(urlFormat(2:end))];
	%urlFormat = '&DELIVERY_FORMAT=CDF';        % default is CDF (3.3) format
    urlDeliveryInterval = '&DELIVERY_INTERVAL=ALL';
	cfaServer = 'http://cfadev.esac.esa.int/'; % default server    
    downloadDirectory = './CFA/';              % local directory where to put downloaded data, default in current directory under 'CFA' subdirectory
    % add
    retrievalType = '&PRODUCT';         % default is to download, not check inventory, hmmmmmm
    noBrowser   = '&NO_BROWSER';        % default is to not download through browser    
    cfaUrl = 'cfa/aio/product-action?'; % used for data retrieval
    selectedFields = 'SELECTED_FIELDS=DATASET_INVENTORY&RESOURCE_CLASS=DATASET_INVENTORY'; % for inventory checks    
    if ~exist('queryFormat','var'); queryFormat = '&RETURN_TYPE=CSV'; end % CFA query format, it is the most readable in Matlab with disp()
    % url encoding: urlencode.m gets ' ' wrong, so uses these instead
    % some places '%' is written directly as '%25'
    space = '%20'; % space-sign: ' '
    and = [space 'AND' space]; % ' AND '
end

%% Check status of downloads if needed
if doLog && checkDownloadStatus,    % check/show status of downloads from .caa file
	disp('=== status of jobs (saved in file .caa) ====');
	if ~isempty(caa),
		for j=1:length(caa), % go through jobs
			disp([num2str(j) '.' caa{j}.status ' ' caa{j}.dataset '-' caa{j}.tintiso]);
		end
	else
		disp('No active downloads');
		if nargout==1, download_status=1; end
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

%% create time interval needed for urls        
if isnumeric(tint) && (size(tint,2)==2), % assume tint is 2 column epoch
	tintiso=irf_time(tint,'tint2iso');
elseif ischar(tint), % tint is in isoformat
	tintiso=tint;
elseif isempty(tint) % will only list products
else
	disp(tint);
	error('caa_download: unknown tint format');
end

if downloadFromCFA % need t1 and t2 instead of t1/t2
    divider=strfind(tintiso,'/');
    t1iso = tintiso(1:divider-1);
    t2iso = tintiso(divider+1:end);
end

%% expand wildcards
if expandWildcards, % expand wildcards
	dataset(strfind(dataset,'?'))='*'; % substitute  ? to * (to have the same convention as in irf_ssub)
	if (any(strfind(dataset,'CIS')) || any(strfind(dataset,'CCODIF')) || any(strfind(dataset,'HIA')))
		dataset(strfind(dataset,'_'))='*'; % substitute  _ to * (to handle CIS products that can have - instead of _)
    end   
end

%% list data if required
if strfind(dataset,'list'),     % list files
	if strcmpi(dataset,'list') && strcmpi(dataset,'listdesc'), % list all files
		filter='*';
	else                        % list only filtered files
		filter=dataset(strfind(dataset,':')+1:end);
	end
	if strfind(dataset,'listdata')
		ttTemp = caa_download(['list:' filter]);
		tint = ttTemp.TimeInterval(1,:);
		ttTemp = caa_download(tint,['list:' filter]);
		iData=find([ttTemp.UserData(:).number]);
		download_status=select(ttTemp,iData);
		return
	end
	if isempty(tint) % work on all datasets
		if any(strfind(dataset,'listdesc')) || any(strfind(dataset,'listgui'))	% get also description
			url_line_list=[caaQuery urlIdentity '&dataset_id=' filter '&dataset_list=1&desc=1'];
			returnTimeTable='listdesc';
		else							% do not get description
			url_line_list=[caaQuery urlIdentity '&dataset_id=' filter '&dataset_list=1'];
			returnTimeTable='list';
		end
    else
        switch downloadFromCFA
            case 0 % CAA
        		url_line_list=[caaInventory urlIdentity '&dataset_id=' filter '&time_range=' tintiso ];
                returnTimeTable='inventory';
            case 1 % CFA
                cfaAction = 'cfa/aio/metadata-action?';
                queryData = ['&QUERY=DATASET.DATASET_ID' space 'like' space '''' cfaQueryDataset ''''];
                queryTime = [and 'DATASET_INVENTORY.START_TIME' space '<=' space '''' t2iso '''',...
                             and 'DATASET_INVENTORY.END_TIME' space '>=' space '''' t1iso ''''];        
                url_line_list = [cfaServer cfaAction selectedFields,...
                         queryData queryTime queryFormat];                 
        end
	end
	disp('Be patient! Contacting CAA...');
	disp(url_line_list);
	caalog=urlread(url_line_list);
	if strfind(dataset,'listgui'), % make gui window with results
		B=regexp(caalog,'(?<dataset>[C][-\w]*)\s+(?<tint>\d\d\d\d-\d\d-\d\d\s\d\d:\d\d:\d\d\s\d\d\d\d-\d\d-\d\d\s\d\d:\d\d:\d\d)\t(?<title>[^\n\t]*)\t(?<description>[^\n\t]*)\n','names');
		list={B.dataset};
		values=cell(numel(B),3);
		values(:,1)={B.tint}';
		values(:,2)={B.title}';
		values(:,3)={B.description}';
		caa_gui_list(list,values)
	else
		if nargout == 1,
			if isempty(returnTimeTable),
				out = textscan(caalog, '%s', 'delimiter', '\n'); % cell array with lines
				download_status = out{1};
			else
				download_status = construct_time_table(caalog,returnTimeTable);
			end
		else
			disp(caalog);
		end
	end
	return;
end
 
%% download data
% create CAA/CFA directory if needed
switch downloadFromCFA
    case 0
        if ~exist('CAA','dir'), mkdir('CAA');end
    case 1
        if ~exist('CFA','dir'), mkdir('CFA');end
end

if checkIfDataAreAtCaa
    switch downloadFromCFA
        case 0 % CAA
            url_line_list=[ caaInventory urlIdentity '&dataset_id=' dataset '&time_range=' tintiso];
            disp(url_line_list);
            disp('Be patient! Contacting CAA to see the list of files...');
            caalist=urlread(url_line_list);
            disp(caalist);
            if ~any(strfind(caalist,'Version')) % there are no CAA datasets available
                disp('There are no CAA data sets available!');
                return;
            end
        case 1 % CFA            
            cfaAction = 'cfa/aio/metadata-action?';
            queryData = ['&QUERY=DATASET.DATASET_ID' space 'like' space '''' cfaQueryDataset ''''];
            queryTime = [and 'DATASET_INVENTORY.START_TIME' space '<=' space '''' t2iso '''',...
                         and 'DATASET_INVENTORY.END_TIME' space '>=' space '''' t1iso ''''];        
            url_line_list = [cfaServer cfaAction selectedFields,...
                         queryData queryTime queryFormat]; 
            disp(url_line_list);
            disp('Be patient! Contacting CFA to see the list of files...');
            caalist=urlread(url_line_list);
            disp(caalist);
            if isempty(caalist) % there are no CFA datasets available
                disp('There are no CFA data sets available!');
                return;
            end
    end  
end
 
switch downloadFromCFA
    case 0 % CAA
        url_line=[caaQuery urlIdentity '&dataset_id=' ...
            dataset '&time_range=' tintiso urlFormat ... 
            urlFileInterval urlNonotify urlSchedule];
        disp('Be patient! Submitting data request to CAA...');
    case 1 % CFA
        url_line = [cfaServer cfaUrl urlIdentity ... 
            '&DATASET_ID=' dataset '&START_DATE=' t1iso '&END_DATE=' t2iso ...
             urlFormat urlDeliveryInterval '&NON_BROWSER' urlNonotify];
        disp('Be patient! Submitting data request to CFA...');
end

disp(url_line);

% Problem with CFA, can only download one file at a time... How to
% implement this? They added this option now.
 
[status,downloadedFile] = get_zip_file(url_line);
if nargout==1, download_status=status;end
if status == 0 && exist(downloadedFile,'file')
	irf_log('fcal','Could not find zip file with data! ');
	fid=fopen(downloadedFile);
	while 1
		tline = fgetl(fid);
		if ~ischar(tline), break, end
		disp(tline)
		if any(strfind(tline,'http:')) && any(strfind(tline,'zip')),
			downloadfile = tline(strfind(tline,'http:'):strfind(tline,'zip')+3);
		end
	end
	fclose(fid);
	delete(downloadedFile);
	
	if exist('downloadfile','var'),
		if doLog
			j=length(caa)+1;
			caa{j}.url=url_line;
			caa{j}.dataset=dataset;
			caa{j}.tintiso=tintiso;
			caa{j}.zip = downloadfile;
			caa{j}.status = 'SUBMITTED';
			caa{j}.timeofrequest = now;
			disp('=====');
			disp('The request has been put in queue');
			disp(['When ready data will be downloaded from: ' downloadfile]);
			disp('To check the status of jobs execute: caa_download');
			caa_log({'Request put in queue: ',url_line,...
				'When ready download from:' downloadfile});
			save -mat .caa caa
		end
		if nargout>=1, download_status=0; end	% 0 if job submitted
	else
		if doLog
			disp('!!!! Did not succeed to download !!!!!');
			caa_log('Did not succeed to download');
		end
		download_status=[];
	end
end
 
%% Functions
function [status,downloadedFile]=get_zip_file(urlLink)
% download zip file, if succeed status=1 and file is unzipped and moved
% to data directory, downloadedFile is set to empty. If there is no zip
% file or file is not zip file, status=0 and downloadedFile is set to
% the downloaded file. 
status = 0; % default
switch downloadFromCFA
    case 0 % CAA
        [downloadedFile,isZipFileReady]=urlwrite(urlLink,tempname);
    case 1 % CFA
        fileName=tempname;
        gzFileName = [fileName '.gz'];
        [gzFileName,isZipFileReady]=urlwrite(urlLink,gzFileName);                
        downloadedFile = gzFileName;
end            		

if isZipFileReady, %
    irf_log('dsrc',['Downloaded: ' urlLink]);
    irf_log('dsrc',['into ->' downloadedFile]);
    caa_log({'Zip file returned for request',urlLink});
    tempDirectory=tempname;
    mkdir(tempDirectory);
    try
        switch downloadFromCFA
            case 0
                filelist=unzip(downloadedFile,tempDirectory);
            case 1
                gunzip(gzFileName);
                filelist=untar(fileName,tempDirectory);                
        end
        if isempty(filelist)
            irf_log('dsrc','Returned zip file is empty');
            caa_log('Zip file empty.');
        else
            move_to_caa_directory(filelist);
        end
        status=1;
        delete(downloadedFile);
        downloadedFile = '';
    catch
        irf_log('fcal','Invalid zip file')
    end
    rmdir(tempDirectory,'s');
else
    irf_log('dsrc',['There is no zip file: ' urlLink]);
end
end
function move_to_caa_directory(filelist)
for jj=1:length(filelist),
    isDataSet = ~any(strfind(filelist{jj},'log'));
    if isDataSet, % dataset files (cdf_convert_summary.log not copied)
        ii=strfind(filelist{jj},filesep);
        dataset=filelist{jj}(ii(end-1)+1:ii(end)-1);
        datasetDirName = [downloadDirectory dataset];
        if ~exist(datasetDirName,'dir'),
            irf_log('dsrc',['Creating directory: ' datasetDirName]); 
            mkdir(datasetDirName);
        elseif overwritePreviousData
            delete([datasetDirName filesep '*']);
        end
        irf_log('dsrc',['file:      ' filelist{jj}]);
        irf_log('dsrc',['moving to directory: ' datasetDirName]);
        movefile(filelist{jj},datasetDirName);
    end
end
end
function paramOut=urlparameter(paramIn)
    if paramIn(1)~= '&'
        paramOut=['&' paramIn];
    else
        paramOut=paramIn;
    end;
end
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
                irf_log('fcal','log file cannot be opened, no log entry');
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
function queryDataset = cfaQueryDataset
    % for wildcards, inventory requests use '%' as wildcard, 
    % while data requests use '*' (something that was not easy to implement)
    if strfind(dataset,'list:')
        queryDataset = dataset(6:end); % assumes 5 first chars are 'list:'
    else
        queryDataset = dataset;
    end
    wildcardIndex = strfind(queryDataset,'*');
    queryDataset(wildcardIndex) = '%';
    queryDataset = urlencode(queryDataset);
    % urlencoding.m works fine with '%', maybe tweak it a little bit so
    % that ' ' becomes '%20' as it should, instead of '+'
end
end
function urlIdentity = get_url_identity(archive)
switch archive
    case 'caa'
        caaUser = datastore('caa','user');
        if isempty(caaUser)
            caaUser = input('Input caa username [default:vaivads]:','s');
            if isempty(caaUser),
                disp('Please register at http://caa.estec.esa.int and later use your username and password.');
                caaUser='vaivads';
            end
            datastore('caa','user',caaUser);
        end
        caaPwd = datastore('caa','pwd');
        if isempty(caaPwd)
            caaPwd = input('Input caa password [default:caa]:','s');
            if isempty(caaPwd), caaPwd='caa';end
            datastore('caa','pwd',caaPwd);
        end
        urlIdentity = ['uname=' caaUser '&pwd=' caaPwd];
    case 'cfa' % just duplicate, but for cfa
        cfaUser = datastore('cfa','user');
        if isempty(cfaUser)
            cfaUser = input('Input caa username [default:avaivads]:','s');
            if isempty(cfaUser),
                disp('Please register at ______? and later use your username and password.');
                cfaUser='avaivads';
            end
            datastore('cfa','user',cfaUser);
        end
        cfaPwd = datastore('cfa','pwd');
        if isempty(cfaPwd)
            cfaPwd = input('Input cfa password [default:!kjUY88lm]:','s');
            if isempty(cfaPwd), cfaPwd='!kjUY88lm';end
            datastore('cfa','pwd',cfaPwd);
        end
        urlIdentity = ['USERNAME=' cfaUser '&PASSWORD=' cfaPwd];
end
end
function TT=construct_time_table(caalog,returnTimeTable)
TT=irf.TimeTable;
switch returnTimeTable
	case 'inventory'
		textLine=regexp(caalog,'(?<dataset>[\w-]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<number>\d+)\s*(?<version>[-\d]+)','names');
		startIndices=regexp(caalog,'(?<dataset>[\w-]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<number>\d+)\s*(?<version>[-\d]+)','start');
		TT.UserData(numel(textLine)).dataset = textLine(end).dataset;
		[TT.UserData(:).dataset]=deal(textLine(:).dataset);
		c=num2cell(str2num(strvcat(textLine(:).number)));
		[TT.UserData(:).number]=deal(c{:});
		c=num2cell(str2num(strvcat(textLine(:).version)));
		[TT.UserData(:).version]=deal(c{:});
	case 'list'
		textLine=regexp(caalog,'(?<dataset>[\w-]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n]*)','names');
		startIndices=regexp(caalog,'(?<dataset>[\w-]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n]*)','start');
		TT.UserData(numel(textLine)).dataset = textLine(end).dataset;
		[TT.UserData(:).dataset]=deal(textLine(:).dataset);
		[TT.UserData(:).title] = deal(textLine(:).title);
	case 'listdesc'
		textLine=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n\t]*)\t(?<description>[^\n]*)','names');
		startIndices=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n\t])*\t(?<description>[^\n]*)','start');
		TT.UserData(numel(textLine)).dataset = textLine(end).dataset;
		[TT.UserData(:).dataset]=deal(textLine(:).dataset);
		[TT.UserData(:).title] = deal(textLine(:).title);
		[TT.UserData(:).description] = deal(textLine(:).description);
	otherwise
		return;
end
tintiso=[vertcat(textLine(:).start) repmat('/',numel(startIndices),1) vertcat(textLine(:).end)];
tint=irf_time(tintiso,'iso2tint');
TT.TimeInterval=tint;
TT.Header = strread(caalog(1:startIndices(1)-1), '%s', 'delimiter', sprintf('\n'));
TT.Comment=cell(numel(TT),1);
TT.Description=cell(numel(TT),1);
end
