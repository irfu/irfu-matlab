function [TS,dobj] = cn_get_ts(filePrefix,dataVar,tt)
% Load mms data file that contains the given time.
%   [ts,dobj] = CN_GET_TS(filePrefix,fileVariable,time);
%   filePrefix - e.g. mms3_edp_brst_l2_scpot
%   fileVariable - e.g. mms3_edp_scpot
%   t - time in TSeries format

isDay = 0;

% Time representation
utc = tt.toUtc;
t.year  = str2double(utc(1:4));
t.month = str2double(utc(6:7));
t.day   = str2double(utc(9:10));
t.hour  = str2double(utc(12:13));
t.min   = str2double(utc(15:16));
t.sec   = str2double(utc(18:end-1));
tdateNum = datenum(t.year,t.month,t.day,t.hour,t.min,t.sec);
disp(['Time: ' utc])

% Directory information
C = strsplit(filePrefix,'_');
varDir = '/data/mms';
for ix=1:length(C), varDir = [varDir filesep C{ix}]; end
dateDir = [utc(1:4) filesep utc(6:7) filesep utc(9:10)];
totDir = [varDir filesep dateDir];


% List files in directory
listingD = dir([varDir filesep dateDir filesep filePrefix '*.cdf']);
if isempty(listingD)
    listingD = dir([varDir filesep dateDir(1:end-3) filesep filePrefix '*.cdf']);
    if isempty(listingD), disp('Could not find any files during the same day!'); TS = []; dobj = []; return;
    else, totDir = totDir(1:end-3); isDay = 1; end
end
disp(['File directory: ' totDir])

% Find file that contains interval start time, t1
nFiles = numel(listingD);
for ii = 1:nFiles, 
    disp(listingD(ii).name)
    if isDay, dateStr = [listingD(ii).name(end-18:end-11) '000001'];
    else dateStr = listingD(ii).name(end-24:end-11);
    end    
    dateNums(ii,1) = datenum(str2double(dateStr(1:4)),str2double(dateStr(5:6)),str2double(dateStr(7:8)),str2double(dateStr(9:10)),str2double(dateStr(11:12)),str2double(dateStr(13:14))); 
end
% plot(1:nFiles,dateNums,'-*',[1 nFiles],tdateNum*[1 1])
fileId = find(tdateNum>dateNums,1,'last');
if isempty(fileId); 
    fileId = 1; 
    allId = fileId;
else
    lastId = fileId;
    allId = find(dateNums==dateNums(lastId));
end
nId = numel(allId);
if nId == 1
    disp('Found 1 file!')
    disp(['File: ' listingD(fileId).name ' (' num2str(listingD(fileId).bytes) ' bytes)'])
else
    disp(['Found ' num2str(nId) ' files!'])
    for ii = 1:nId        
        disp(['File ' num2str(ii) ': ' listingD(allId(ii)).name ' (' num2str(listingD(allId(ii)).bytes) ' bytes)'])
    end
    specId = input(['Choose file #(' num2str(1) '-' num2str(nId) '): ']);
    fileId = allId(specId);
    disp(['File: ' listingD(fileId).name ' (' num2str(listingD(fileId).bytes) ' bytes)'])
end

%% Load the data
dobj=dataobj([totDir filesep listingD(fileId).name]);
if isempty(dataVar)
    TS = [];
else    
    if ischar(dataVar)
        dataVar = {dataVar};
    end
    for ii = 1:numel(dataVar)
        try   
            ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));
        catch    
            ts = get_variable(dobj,dataVar{ii});
        end    
        TS(ii) = ts;
    end
end