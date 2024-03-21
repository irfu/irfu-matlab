function [maneuvers, timelineXML, eclipse] = mms_maneuvers( Tint, scIdStr )
% MMS_MANEUVERS get maneuver and eclipse information
%	[maneuvers, timelineXML, eclipse] = MMS_MANEUVERS(Tint, [scIdStr]) reads the
%   appropriate FDOA timeline files and extract all maneuvers and eclipses
%   that took place during the time interval "Tint".
%   If "scIdStr" is specified, then it only returns the manuevers / eclipses that
%   occured during the time interval "Tint" on the spacecraft "scIdStr".
%
%   Note: This function relies on MMS_READ_TIMELINE and therefor indirectly
%   on XPath, should FDOA change their file structure this function may
%   require changes.
%
%   Input:
%     Tint    = Time interval to look for manuevers, output from irf.tint().
%     scIdstr = string of which spacecraft to be extracted, values '1',
%               '2', '3' or '4'. Or a combination such as '1234'.
%   Output:
%     manuevers    = struct containing
%       .mms1      = cell of maneuvers on MMS1, in irf.tint time interval.
%       .mms2      = corresponding on MMS2.
%       .mms3      = corresponding on MMS3.
%       .mms4      = corresponding on MMS4.
%     timelineXML  = list of timeline xml files containing the maneuvers.
%                    Note: Only files containing maneuvers, not a complete
%                    list of all xml files read if they do not contain any
%                    relevant manuevers for the requested interval.
%   Optional output:
%     eclipse     = struct similar to maneuvers but intervals of eclipses
%                   instead of maneuvers.
%
%   Example:
%
%   % Setup time interval
%   Tint = irf.tint('2015-03-14T22:50:00Z/2016-04-01T23:00:00Z');
%	maneuvers = MMS_MANEUVERS(Tint, '1');
%    returns a struct maneuvers containing
%       .mms1 = Cell with time intervals of all maneuvers for MMS1,
%               (in irf.tint format).
%       .mms2 = Empty cell (scIdStr was specified to only look for MMS1)
%       .mms3 = Empty cell
%       .mms4 = Empty cell
%
% 	See also MMS_READ_TIMELINE.

% Verify input
narginchk(1,2);
if(~isa(Tint,'GenericTimeArray') || ~all(size(Tint)==[2,1])), error('Unexpected Tint'); end
if(nargin<2), scIdStr = '1234'; end
% Ensure output
maneuvers.mms1={}; maneuvers.mms2={}; maneuvers.mms3={}; maneuvers.mms4={};
if(nargout>=3)
  lookForEclipse=true;
  eclipse.mms1={}; eclipse.mms2={}; eclipse.mms3={}; eclipse.mms4={};
else
  lookForEclipse=false;
end

% List timeline files (can cover up to over 100 days...) The following is a
% quick fix, should be integrated into mms db list_files.
% /data/mms/ancillary/mms/timeline/mms_timeline_yyyyDOY_yyyyDOY_vXX.xml
dataRootPath = getenv('DATA_PATH_ROOT');
if(isempty(dataRootPath) || ~exist(dataRootPath,'dir'))
  dataRootPath='/data/mms'; % Fallback for IRFU location of data
end
dataPath = [dataRootPath, filesep, 'ancillary', filesep, 'mms', ...
  filesep, 'timeline', filesep];
% List all timeline xml files.
list = dir([dataPath,'mms_timeline_*.xml']);
for ii = length(list):-1:2
  % If the filename of the next entry match the present file (excl. version
  % number) remove the superseeded file from list.
  if(strcmpi(list(ii-1).name(1:end-6), list(ii).name(1:end-6)))
    irf.log('debug',['Removing superseeded file ',list(ii-1).name,' from list.']);
    list(ii-1)=[];
  end
end

% Use a copy of Tint, so that it can be reduced as reading progress.
Tint2 = Tint;

% Read list from End to beginning, if Time of file is within interval Tint2,
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
    irf.log('notice',['Processing file: ',list(ii).name]);
    try
      if lookForEclipse
        [maneuv, fileInterval, eclip] = mms_read_timeline([dataPath, list(ii).name], scIdStr);
      else
        [maneuv, fileInterval] = mms_read_timeline([dataPath, list(ii).name], scIdStr);
      end
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
      if(isfield(maneuv, ['mms',num2str(kk)]))
        tmp = maneuv.(['mms',num2str(kk)]);
        for ll=length(tmp):-1:1
          % Check if it is outside of Tint2
          if( tmp{ll}.start.ttns >= Tint2.stop.ttns || ...
              tmp{ll}.stop.ttns < Tint2.start.ttns )
            irf.log('debug',['Skipping maneuver outside of interval ', ...
              Tint2.start.toUtc(1), '/', Tint2.stop.toUtc(1)]);
            continue
          else
            % Append
            irf.log('debug', ['Maneuver on MMS', num2str(kk), ' ', ...
              tmp{ll}.start.toUtc(1), '/', tmp{ll}.stop.toUtc(1), ...
              ' inside of requested interval.']);
            if(~isempty(maneuvers.(['mms',num2str(kk)])))
              tmp2 = maneuvers.(['mms',num2str(kk)]);
              if(tmp2{1}.start.ttns ~= tmp{ll}.start.ttns)
                maneuvers.(['mms',num2str(kk)]) = [tmp(ll); tmp2];
                tmpTimeline = [tmpTimeline; list(ii).name]; %#ok<AGROW>
              else
                irf.log('debug', 'Maneuver already included.');
              end
            else
              maneuvers.(['mms',num2str(kk)]) = tmp(ll);
              tmpTimeline = list(ii).name;
            end
          end
        end % for ll
      end % isfield
      if(lookForEclipse)
        % Eclipse
        if(isfield(eclip, ['mms',num2str(kk)]))
          tmp = eclip.(['mms',num2str(kk)]);
          for ll=length(tmp):-1:1
            % Check if it is outside of Tint2
            if( tmp{ll}.start.ttns >= Tint2.stop.ttns || ...
                tmp{ll}.stop.ttns < Tint2.start.ttns )
              irf.log('debug',['Skipping eclipse outside of interval ', ...
                Tint2.start.toUtc(1), '/', Tint2.stop.toUtc(1)]);
              continue
            else
              % Append
              irf.log('debug', ['Eclipse on MMS', num2str(kk), ' ', ...
                tmp{ll}.start.toUtc(1), '/', tmp{ll}.stop.toUtc(1), ...
                ' inside of requested interval.']);
              if(~isempty(eclipse.(['mms',num2str(kk)])))
                tmp2 = eclipse.(['mms',num2str(kk)]);
                if(tmp2{1}.start.ttns ~= tmp{ll}.start.ttns)
                  eclipse.(['mms',num2str(kk)]) = [tmp(ll); tmp2];
                  tmpTimeline = [tmpTimeline; list(ii).name]; %#ok<AGROW>
                else
                  irf.log('debug', 'Eclipse already included.');
                end
              else
                eclipse.(['mms',num2str(kk)]) = tmp(ll);
                tmpTimeline = list(ii).name;
              end
            end
          end % for ll
        end % isfield
      end % Eclipse
    end % for kk
    % Replace Tint2 interval to look for maneuvers in next file to cover
    % times up to the start of this previous file (ie. events that have
    % been changed or cancelled should not be included from older files).
    tEnd = max(min(fileInterval.start.ttns,Tint2.stop.ttns), Tint2.start.ttns);
    Tint2 = irf.tint(Tint2.start, EpochTT(tEnd));
    if(Tint2.start.ttns - Tint2.stop.ttns == int64(0))
      % We have processed all of the latest timeline files that cover our
      % original requested time interval Tint.
      % Improve?: Possibly limit resulting maneuvers that cover longer time
      % intervals than the specified original Tint.

      % if maneuvers found and tmpTimeline exist, return only the unique
      % list of files
      if(exist('tmpTimeline','var') && ~isempty(tmpTimeline))
        timelineXML = unique(tmpTimeline, 'rows');
      else
        timelineXML = [];
      end
      return
    end
  end
end

end
