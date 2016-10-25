function varargout = get_filepath(filePrefix,time,varargin)
% MMS.get_filepath Get path to mms data file.
%   Simple finding routine that only checks that the file name should 
%   contain the given time interval. You need to verify that the file 
%   indeed contains the expected time interval manually.
%
%   Before running, use mms.db_init('local_file_db','path_to_database');
%     for example: mms.db_init('local_file_db','/Volumes/ExternalDrive/data');
%
%   filepath_and_filename = MMS.GET_FILEPATH(filePrefix,time);
%       filePrefix - e.g. mms3_edp_brst_l2_scpot
%       time - time or time interval in EpochTT format
%
%   [filepath,filename] = MMS.GET_FILEPATH(filePrefix,time);
%       filepath = '/Volumes/ExternalDrive/data/mms1/fpi/brst/l2/des-dist/2015/10/16/'
%       filename = 'mms1_fpi_brst_l2_des-dist_20151016102944_v2.1.0.cdf'
%
%   MMS.GET_FILEPATH(filePrefix,time,'list'); 
%       - displays list of files in the folder corresponding to the time given
%         
%   MMS.GET_FILEPATH(filePrefix,time,'noroot'); 
%       - exclude the root of the filepath e.g. '/Volumes/ExternalDrive/data/'
%         so that only relative path is returned, e.g.: '/mms1/fpi/brst/l2/des-dist/2015/10/16/'

isDay = 0;
doList = 0;
noRoot = 0;

% Check for additional options
have_options = 0;
if numel(varargin) > 0, have_options = 1; args = varargin; end
while have_options
  switch(lower(args{1}))   
    case 'list'
      doList = 1;
      args = args(2:end); 
    case 'noroot'
      noRoot = 1;
      args = args(2:end); 
  end 
  if isempty(args), break, end    
end

%if nargin == 3 && strcmp(varargin{1},'list'); doList = 1; end

% Time representation
utc = time(1).toUtc;
t.year  = str2double(utc(1:4));
t.month = str2double(utc(6:7));
t.day   = str2double(utc(9:10));
t.hour  = str2double(utc(12:13));
t.min   = str2double(utc(15:16));
t.sec   = str2double(utc(18:end-1));
tdateNum = datenum(t.year,t.month,t.day,t.hour,t.min,t.sec);
disp(['Time: ' utc])

% Directory information
db_info = datastore('mms_db'); 
%dataDir = '/data/mms';
dataDir = db_info.local_file_db_root;
varDir = strjoin(strsplit(filePrefix,'_'),filesep);
dateDir = [utc(1:4) filesep utc(6:7) filesep utc(9:10)];
totDir = [dataDir filesep varDir filesep dateDir filesep];

% List files in directory
listingD = dir([totDir filesep filePrefix '*.cdf']);
if isempty(listingD) % Files are in month folder
    totDir = totDir(1:end-3); isDay = 1;
    listingD = dir([totDir filesep filePrefix '*.cdf']);
    if isempty(listingD) % Also month folder is empty
        disp(['Could not find any ' filePrefix ' files during the same day!']); 
        TS = []; dobj = []; return;
    end
end
disp(['  Directory: ' totDir])
nFiles = numel(listingD);
if doList % Only list the files
    for ii = 1:nFiles 
        disp(listingD(ii).name)        
    end
    TS = []; dobj = [];
    return;
end

% Find file that contains interval start time, t1
for ii = 1:nFiles     
    splitName = strsplit(listingD(ii).name,'_');
    dateStr = splitName{end-1};
    if isDay, dateStr = [dateStr '000001']; end
    dateNums(ii,1) = datenum(str2double(dateStr(1:4)),str2double(dateStr(5:6)),str2double(dateStr(7:8)),str2double(dateStr(9:10)),str2double(dateStr(11:12)),str2double(dateStr(13:14))); 
end
% plot(1:nFiles,dateNums,'-*',[1 nFiles],tdateNum*[1 1])
fileId = find(tdateNum>dateNums,1,'last');
if isempty(fileId) 
    fileId = 1; 
    allId = fileId;
else
    lastId = fileId;
    allId = find(dateNums==dateNums(lastId));
end
nId = numel(allId);


if nId > 1            
    %specId = input(['Choose file #(' num2str(1) '-' num2str(nId) '): ']);
    %fileId = allId(specId);
    fileId = allId(end); % just load highest version
    %disp('Loading highest version.')    
end
disp(['  File: ' listingD(fileId).name ' (' num2str(listingD(fileId).bytes) ' bytes)'])

if noRoot
  filepath = [varDir filesep dateDir filesep];
else
  filepath = totDir;
end
filename = listingD(fileId).name;

if nargout == 0
  varargout = {[filepath  filename]};
elseif nargout == 1
  varargout = {[filepath  filename]};
elseif nargout == 2
  varargout = {filepath,filename};
end
