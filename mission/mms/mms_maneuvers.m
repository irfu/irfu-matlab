function maneuvers = mms_maneuvers( Tint, scIdStr )
% Read Timeline and return maneuvers in the interval Tint
%
% TESTING: DO NOT USE!

narginchk(1,2);
maneuvers.mms1={}; maneuvers.mms2={}; maneuvers.mms3={}; maneuvers.mms4={};

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
  if( ( (fileStart - 86400*10^9) > Tint2.stop.ttns ) || ...
      ( (fileStop + 86400*10^9) < Tint2.start.ttns ) )
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
    if( fileInterval.start.ttns > Tint2.stop.ttns || ...
        fileInterval.stop.ttns < Tint2.start.ttns )
      % Skip file as it did not cover Tint time inteval despite its name
      % suggested it would..
      irf.log('notice',['Skipping file : ',list(ii).name]);
      continue
    end
    % Locate maneuv within Tint2 (contains time interval from original
    % Tint.start to the previous fileInterval.start) and append these to
    % maneuvers (to be returned).
    for kk=1:4
      tmp = maneuv.(['mms',num2str(kk)]);
      for ll=length(tmp):-1:1
        % Check if it is outside of Tint2
        if( tmp{ll}.start.ttns > Tint2.stop.ttns || ...
            tmp{ll}.stop.ttns < Tint2.start.ttns )
          continue
        else
          % Append
          if(~isempty(maneuvers.(['mms',num2str(kk)])))
            tmp2 = maneuvers.(['mms',num2str(kk)]);
            maneuvers.(['mms',num2str(kk)]) = [tmp(ll); tmp2];
          else
            maneuvers.(['mms',num2str(kk)]) = {tmp{ll}};
          end
        end
      end
    end
    % Replace Tint2 interval to look for maneuvers in next file to cover
    % times up to the start of this previous file (ie. events that have
    % been changed or cancelled should not be included from older files).
    tEnd = max(min(fileInterval.start.ttns,Tint2.stop.ttns), Tint2.start.ttns);
    Tint2 = irf.tint(Tint2.start, EpochTT(tEnd));
    if(Tint2.start.ttns - Tint2.stop.ttns == int64(0)),
      % We have processed all of the latest timeline files that cover our
      % original requested time interval Tint.
      return
    end
  end
end

end