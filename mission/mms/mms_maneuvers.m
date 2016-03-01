function maneuvers = mms_maneuvers( Tint, scIdStr )
% Read Timeline and return maneuvers in the interval Tint
%
% TESTING: DO NOT USE!

narginchk(1,2);
maneuvers = struct('mms1',[],'mms2',[],'mms3',[],'mms4',[]);

if(~isa(Tint,'GenericTimeArray') || ~all(size(Tint)==[2,1])), error('Unexpected Tint'); end
if(nargin<2), scIdStr = '1234'; end

% List timeline files (can cover up to over 100 days...) The following is a
% quick fix, should be integrated into mms db list_files.
% /data/mms/ancillary/mms/timeline/mms_timeline_yyyyDOY_yyyyDOY_vXX.xml
%% FIXME: "mounted/" should be removed when put in GIT.
if(strcmp(getComputerName, 'thonilaptop'))
dataPath = [filesep, 'data', filesep, 'mms', filesep, 'mounted', ...
  filesep, 'ancillary', filesep, 'mms', filesep, 'timeline', filesep];
else
  error('DO NOT USE THIS FUNCTION YET!');
end

% List all timeline files and remove superseeded versions.
list = dir([dataPath,'mms_timeline_*.xml']);
for ii = length(list):-1:2
  % If next file match the current file, remove next file from list.
  if(strcmpi(list(ii-1).name(1:end-6), list(ii).name(1:end-6)))
    list(ii-1)=[];
  end
end

%% 

Tint2 = Tint;

% Read list from End to beginning, if Time of file is within interval Tint,
% and for each file only include times before the creation time of the
% previous file. (Avoid Maneuvers that might have been changed in time or
% cancelled).
for ii = length(list):-1:1
  nameInfo = sscanf(list(ii).name, 'mms_timeline_%4u%3u_%4u%3u_v%u.xml');
  % nameInfo(1) - Start YYYY, nameInfo(2) - Start DOY
  % nameInfo(3) - Stop  YYYY, nameInfo(4) - Stop  DOY, nameInfo(5)-Verison
  fileStart = irf_time([nameInfo(1), nameInfo(2)],'doy>ttns');
  fileStop = irf_time([nameInfo(3), nameInfo(4)], 'doy>ttns');
  if( ( (fileStart - 86400*10^9) > Tint.stop.ttns ) || ...
      ( (fileStop + 86400*10^9) < Tint.start.ttns ) )
    % Remove file from list.
    list(ii) = [];
  else
    % Load list(ii), check which if any maneuvers match our requested Tint.
    % In the next list(ii) add only those maneuvers that match our reqested
    % Tint and that are before the creation date of the previous list(ii)
    % file. This way any maneuvers that have been cancelled or moved in time
    % should not be included in our final list.
    disp(list(ii).name);
    try
      [maneuv, fileInterval] = mms_read_timeline([dataPath, list(ii).name], scIdStr);
    catch ME
      errStr = ['Failed to read file: ', dataPath,list(ii).name, ' with error message: ',...
        ME.message, ' Trying the next file.'];
      irf.log('warning', errStr);
      continue
    end
    
    %% FIXME combine "maneuv" (just read) and "maneuvers" (read in previous file)
    keyboard
    
    % Replace Tint interval to look for events to cover times up to this
    % file for the following files (ie. events that have been
    % changed/cancelled should not be loaded from older files).
    tEnd = max(min(fileInterval.start.ttns,Tint2.stop.ttns), Tint.start.ttns);
    Tint2 = irf.tint(Tint.start, EpochTT(tEnd));
  end
end

end