function [download_status,downloadfile]=caa_download(tint,dataset,varargin)
% CAA_DOWNLOAD Download CAA data in CDF format
%       CAA_DOWNLOAD - check the status of jobs in current directory
%
%       CAA_DOWNLOAD('list')    - list all datasets and their available times
%       TT=CAA_DOWNLOAD('list')- return time table with datasets and available times
%       CAA_DOWNLOAD('listdesc')- same with dataset description
%       CAA_DOWNLOAD('listgui') - same presenting output in separate window
%       CAA_DOWNLOAD('list:dataset')- list:/listdesc:/listgui:  filter datasets 'dataset'
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
%   'file_interval' - see command line manual http://goo.gl/VkkoI, default 'file_interval=72hours'
%   'nowildcard'	- download the dataset without any expansion in the name and not checking if data are there
%   'overwrite'		- overwrite files in directory (to keep single cdf file) NEEDS IMPLEMENTATION
%   'schedule'		- schedule the download, (returns zip file link)
%						check the readiness by executing CAA_DOWNLOAD from the same direcotry
%   'nolog'			- do not log into .caa file (good for batch processing)
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
%   'test'        - use caa test server instead
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%% Defaults
checkDownloadsStatus	= false;
doLog					= true; % log into .caa file
flag_wildcard     =1;                     % default is to use wildcard
flag_check_if_there_is_data=1;            % check if there are any at caa
urlNonotify='&nonotify=1';                % default is not notify by email
urlFileInterval='&file_interval=72hours'; % default time interval of returned files
urlSchedule='';                           % default do not have schedule option
urlFormat='&format=cdf';                  % default is CDF (3.3) format
caaServer='http://caa.estec.esa.int/'; % default server
urlIdentity='?uname=vaivads&pwd=caa';     % default identity
%urlInventory='';                          % default no inventory output (currently use different www link)
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
if nargin==0, checkDownloadsStatus=true; end
if nargin>2, % cehck for additional flags
	for iFlag=1:numel(varargin)
		flag=varargin{iFlag};
		if strcmpi(flag,'test'),  % use test server
			caaServer='http://caa5.estec.esa.int/caa_query/';
			urlNonotify='';           % notify also by email
		elseif strcmpi(flag,'nowildcard'),
			flag_wildcard=0;
			flag_check_if_there_is_data=0;
			urlNonotify='&nonotify=1';
		elseif any(strfind(flag,'file_interval'))
			urlFileInterval = urlparameter(flag);
		elseif any(strcmpi('schedule',flag))
			urlSchedule = '&schedule=1';
		elseif any(strcmpi('nolog',flag))
			doLog = false;
		elseif any(strcmpi('inventory',flag))
			urlSchedule = '&inventory=1';
		else
			irf_log('fcal',['Flag ''' flag ''' not recognized']);
		end
	end
end
if nargin>=1, % check if fist argument is not caa zip file link
	if ischar(tint)
		if nargin>1 && ischar(dataset) && strcmpi(dataset,'nolog')
			doLog=false;
		end
		if regexp(tint,'\.zip') % download data zip file
			if doLog
				j=numel(caa)+1;
				caa{j}.url='*';
				caa{j}.dataset='*';
				caa{j}.tintiso='*';
				caa{j}.zip=tint;
				caa{j}.status='submitted';
				caa{j}.timeofrequest=now;
				checkDownloadsStatus=true;
			else
				temp_file=tempname;
				zipFileLink=tint;
				isJobFinished=get_zip_file(zipFileLink,temp_file);
				if isJobFinished, %
					download_status=1;
					return;
				else
					irf_log('dsrc','Job still not finished');
					download_status=0;
					return;
				end
			end
		else	% list or inventory ingested data, no time interval specified
			dataset=tint;
			tint=[];
		end
	elseif ~isnumeric(tint)
		help caa_download;return;
	end
end
if nargout>0, checkDownloadsStatus=false; end
caaQuery=[caaServer 'caa_query/'];
caaInventory=[caaServer 'cgi-bin/inventory.cgi/'];
%% Check status of downloads if needed
if doLog && checkDownloadsStatus,    % check/show status of downloads from .caa file
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
	j_remove_jobs=zeros(1,length(caa));
	j_finished_jobs=zeros(1,length(caa));
	for j=1:length(caa), % go through jobs
		if strcmpi(caa{j}.status,'downloaded') || strcmpi(caa{j}.status,'finnished') || strcmpi(caa{j}.status,'finished') % 'finnished shoudl be removed after some time % do nothing
			j_finished_jobs(j)=1;
		elseif strcmpi(caa{j}.status,'submitted'),
			disp(['=== Checking status of job nr: ' num2str(j) '==='])
			temp_file=tempname;
			isJobFinished=get_zip_file(caa{j}.zip,temp_file);
			if isJobFinished, %
				caa{j}.status='FINISHED';
				save -mat .caa caa; % changes in caa saved
			else % job not finished
				disp(['STILL WAITING TO FINISH, submitted ' num2str((now-caa{j}.timeofrequest)*24*60,3) 'min ago.']);
				if now-caa{j}.timeofrequest>1, % waiting more than 1 day
					y=input('Waiting more than 24h. Delete from list? y/n :','s');
					if strcmpi(y,'y'),
						j_remove_jobs(j)=1;
					end
				end
			end
		else
			disp('ERROR: Unknown status!')
			return
		end
	end
	if sum(j_finished_jobs)>5, % ask for cleanup
		y=input('Shall I remove FINISHED from the list? y/n :','s');
		if strcmpi(y,'y'),
			j_remove_jobs=j_remove_jobs | j_finished_jobs;
		end
	end
	caa(j_remove_jobs==1)=[];
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

%% expand wildcards
if flag_wildcard, % expand wildcards
	dataset(strfind(dataset,'?'))='*'; % substitute  ? to * (to have the same convention as in irf_ssub)
	if (any(strfind(dataset,'CIS')) || any(strfind(dataset,'CCODIF')) || any(strfind(dataset,'HIA')))
		dataset(strfind(dataset,'_'))='*'; % substitute  _ to * (to handle CIS products that can have - instead of _)
	end
end

%% list data if required
if strfind(dataset,'list'),     % list  files
	if strcmpi(dataset,'list') && strcmpi(dataset,'listdesc'), % list all files
		filter='*';
	else                        % list only filtered files
		filter=dataset(strfind(dataset,':')+1:end);
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
		url_line_list=[caaInventory urlIdentity '&dataset_id=' filter '&time_range=' tintiso ];
		returnTimeTable='inventory';
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
% create CAA directory if needed
if ~exist('CAA','dir'), mkdir('CAA');end

if flag_check_if_there_is_data
	%    url_line_list=['http://caa.estec.esa.int/caa_query/?uname=vaivads&pwd=caa&dataset_id=' ...
	%        dataset '&time_range=' tintiso '&format=cdf&list=1'];
	%    url_line_list=['http://caa.estec.esa.int/caa_test_query/?uname=vaivads&pwd=caa&dataset_id=' ...
	%        dataset '&time_range=' tintiso '&format=cdf&list=1'];
	url_line_list=['http://caa.estec.esa.int/cgi-bin/inventory.cgi/?uname=vaivads&pwd=caa&dataset_id=' dataset '&time_range=' tintiso];
	disp(url_line_list);
	disp('Be patient! Contacting CAA to see the list of files...');
	caalist=urlread(url_line_list);
	disp(caalist);
	if ~any(strfind(caalist,'Version')),% there are no CAA datasets available
		disp('There are no CAA data sets available!');
		return;
	end
end

url_line=[caaQuery urlIdentity '&dataset_id=' ...
	dataset '&time_range=' tintiso urlFormat urlFileInterval urlNonotify urlSchedule];

disp('Be patient! Submitting data request to CAA...');
disp(url_line);

try
	temp_file=tempname;
	get_zip_file(url_line,temp_file);
	if nargout==1, download_status=1;end
catch
	irf_log('fcal','Could not find zip file with data! ');
	fid=fopen(temp_file);
	while 1
		tline = fgetl(fid);
		if ~ischar(tline), break, end
		disp(tline)
		if any(strfind(tline,'http:')) && any(strfind(tline,'zip')),
			downloadfile = tline(strfind(tline,'http:'):strfind(tline,'zip')+3);
		end
	end
	fclose(fid);
	delete(temp_file);
	
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

	function status=get_zip_file(urlLink,temp_file)
		[f,isZipFileReady]=urlwrite(urlLink,temp_file);
		if isZipFileReady, %
			irf_log('dsrc',['Downloaded: ' urlLink]);
			irf_log('dsrc',['into ->' temp_file]);
			caa_log({'Zip file returned for request',urlLink});
			tempDirectory=tempname;
			filelist=unzip(temp_file,tempDirectory);
			if isempty(filelist)
				irf_log('dsrc','Returned zip file is empty');
				caa_log('Zip file empty.');
			else
				move_to_caa_directory(filelist);
			end
			rmdir(tempDirectory,'s');
			delete(f);
			status=1;
		else
			irf_log('dsrc',['There is no zip file: ' urlLink]);
			status=0;
		end			
	end
	function move_to_caa_directory(filelist)
		for jj=1:length(filelist),
			ii=strfind(filelist{jj},filesep);
			isDataSet = ~any(strfind(filelist{jj},'log'));
			if isDataSet, % dataset files (cdf_convert_summary.log not copied)
				dataset=filelist{jj}(ii(end-1)+1:ii(end)-1);
				disp(['Data set: ' dataset '--> CAA/']);
				if ~exist(['CAA/' dataset],'dir'), mkdir(['CAA/' dataset]);end
				movefile(filelist{jj},['CAA/' dataset]);
			end
		end
		%disp(['REMOVING DATA DIRECTORIES & FILES: ' filelist{jj}(1:ii(1)) ',delme.zip']);
		%rmdir(filelist{jj}(1:ii(1)),'s');
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
end
function TT=construct_time_table(caalog,returnTimeTable)
TT=irf.TimeTable;
switch returnTimeTable
	case 'inventory'
		textLine=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<number>\d+)\s*(?<version>[-\d]+)','names');
		startIndices=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<number>\d+)\s*(?<version>[-\d]+)','start');
		TT.UserData(numel(textLine)).dataset = textLine(end).dataset;
		[TT.UserData(:).dataset]=deal(textLine(:).dataset);
		c=num2cell(str2num(strvcat(textLine(:).number))); 
		[TT.UserData(:).number]=deal(c{:});
		c=num2cell(str2num(strvcat(textLine(:).version))); 
		[TT.UserData(:).version]=deal(c{:});
	case 'list'
		textLine=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n]*)','names');
		startIndices=regexp(caalog,'(?<dataset>[\w]*)\s+(?<start>[\d-]{10}\s[\d:]+)\s*(?<end>[\d-]+\s[\d:]+)\s*(?<title>[^\n]*)','start');
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
TT.Header=caalog(1:startIndices(1)-1);
TT.Comment=cell(numel(TT),1);
TT.Description=cell(numel(TT),1);
end
