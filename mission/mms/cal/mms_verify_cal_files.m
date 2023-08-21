% Function to verify the latest mms cal file is correctly following the
% expected format

function calOk = mms_verify_cal_files()

calOk = true;
calPath = [irf('path'), filesep, 'mission', filesep, 'mms', filesep, 'cal'];
calTypes = {'dsl', 'scpot', 'regions'};
for calId = 1:length(calTypes)
  calStr = calTypes{calId};
  for scId=1:4
    try
      calFiles = [calPath, filesep, 'mms', num2str(scId), '_edp_sdp_', ...
        calStr, '_*v*.txt' ];
      list = dir(calFiles);
      if(~isempty(list))
        list = list(end); % Check the last version, ie. end of list.
        verStr = regexp(list.name, ['mms[1-4]_edp_sdp_', calStr, ...
          '_\d{8,8}_v(?<VERSION>\d{1,}.\d{1,}.\d{1,})'],'tokens');
        if is_version_geq(verStr{1}{1}, '1.0.0') % Version check
          % NEW format of Calibration files
          % As discussed during meeting at KTH, 2017/10/?? the new calibration
          % files have two extra columns providing the time margins to be
          % applied before and after each timestamp.
          %         __________
          %        /
          %      ./
          %_____/
          %
          %      | <-- Times given by first column, but the value does not come
          %      into full effect until after the margin of the last column.
          %      But it starts to come into effect before the margin of the
          %      second to last column. If margins are zero, then the old value
          %      is used until the timestamp in the first column and after this
          %      timestamp the new offsets are used. Otherwise a linear
          %      interpolation is applied in the intermediate interval, before
          %      which the previous offset is used and correspondingly after
          %      the margin the new offset is fully used.
          fID = fopen([calPath, filesep, list.name]);
          C = textscan(fID, '%s\t%f\t%f\t%f\t%f', 'HeaderLines', 1);
          fclose(fID);
          %% Verify the new verison "dt" margins
          [ind, ~] = find(C{4} > 0); % all "-dt" should be negative or zero
          if any(ind)
            warning('Found %i rows with non-negative time margin "-dt", %s', length(ind), list.name);
            warning(['First (max 5) rows with problmes were: ' sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
            calOk = false;
          end
          [ind, ~] = find(C{5} < 0); % all "+dt" should be positive or zero
          if any(ind)
            warning('Found %i rows with non-positive time margin "+dt", %s', length(ind), list.name);
            warning(['First (max 5) rows with problmes were: ' sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
            calOk = false;
          end


        else
          % OLD format of Calibration files
          fID = fopen([calPath, filesep, list.name]);
          switch calStr
            case 'regions'
              fmt = '%s\t%f';
            case {'dsl', 'scpot'}
              fmt = '%s\t%f\t%f';
            otherwise
              error('Not yet defined type and format');
          end
          C = textscan(fID, fmt, 'HeaderLines', 1);
          fclose(fID);
        end
        % Interpret the time strings and convert it to int64 (tt2000).
        offTime = EpochTT(cell2mat(C{1}));
        % Offset start time (based on ROI).
        time1 = offTime.ttns;

        %% Verify Time is monotone increasing
        [ind, ~] = find(diff(time1)<=0);
        if any(ind)
          warning('Found %i rows with time not being monotonically increasing, %s', length(ind), list.name);
          warning(['First (max 5) rows with problems were: ', sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
          calOk = false;
        end

        %% Verify Time (with "-/+dt") is also monotone increasing
        if is_version_geq(verStr{1}{1}, '1.0.0')
          t_start = time1 + int64(C{4}*10^9); % t -dt
          t_start = t_start(2:end); % Skip start point
          t_stop = time1 + int64(C{5}*10^9); % t +dt
          t_stop = t_stop(1:end-1); % Skip end point
          t_12 =  reshape([t_stop'; t_start'], [], 1);
          [ind, ~] = find(diff(t_12)<=0);
          if any(ind)
            warning('Found %i rows with time not being monotone increasing when considering "+/-dt", %s', length(ind), list.name);
            warning(['First (max 5) rows with problems were: ', sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
            calOk = false;
          end
        end

        %% Verify no NaN values in the offsets
        if(~strcmp(calStr, 'regions'))
          dataOff = [C{2}, C{3}];
        else
          % regions files does not have a third column (as of v0.0.0)
          dataOff = C{2};
        end
        [ind, ~] = find(isnan(dataOff));
        if any(ind)
          ind = unique(sort(ind));
          warning('Found %i rows with NaN in offset file, %s', length(ind), list.name);
          warning(['First (max 5) rows with problems were: ', sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
          calOk = false;
        end

        %% Verify no Inf values in the offsets
        [ind, ~] = find(isinf(dataOff));
        if any(ind)
          ind = unique(sort(ind));
          warning('Found %i rows with Inf in offset file, %s', length(ind), list.name);
          warning(['First (max 5) rows with problems were: ', sprintf('%i, ', 1+ind(1:min(5, length(ind))))]);
          calOk = false;
        end

        %% Verify end value (of region files) is not S/W
        % as we should not have it be active for all time (until next
        % region file is ready for production)
        if strcmp(calStr, 'regions')
          if dataOff(end) == 1
            warning('Region file %s appears to have solar wind active for the last record.', list.name);
            calOk = false;
          end
        end

      else
        warning('Failed to locate %s calibration file for scId: %i', calStr, scId);
        calOk = false;
      end
    catch ME
      if exist('list','var') && ~isempty(list)
        warning('Failed to process file: %s', list.name);
      elseif exist('calFiles', 'var')
        warning('Failed to process file matching: %s', calFiles);
      end
      warning(ME.message);
      calOk = false;
    end
  end
end

end