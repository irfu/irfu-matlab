function indSW = mms_sdp_swwake_enabled_time(tInp, scId)
% Returns index "indSW" of when mms S/W wake removal code should run, and
% when it should not.

narginchk(2,2)
if isa(tInp,'GenericTimeArray'), tInp = EpochTT(tInp).ttns; end
global ENVIR

% NEW CODE, READ FILES PROVIDED BY DANIEL's code (DOES NOT YET EXIST in
% $CAL_PATH_ROOT)...
if isempty(ENVIR), ENVIR = mms_sdc_sdp_init; end
calStr = 'regions';
calPath = [ENVIR.CAL_PATH_ROOT, filesep 'mms', num2str(scId), filesep, ...
  'edp', filesep, 'sdp', filesep, calStr];
calFiles = ['mms', num2str(scId), '_edp_sdp_', calStr, '_*v*.txt' ];
list = dir([calPath, filesep, calFiles]);
if(~isempty(list))
  % Found at least one cal file. Use the last version, ie. end of list.
  list = list(end);
  [~, off.calFile, ~] = fileparts(list.name);
  irf.log('notice', ['Reading S/W times from: ',off.calFile]);
  fmt = '%s\t%f';
  fileID = fopen([calPath, filesep, list.name]);
  C = textscan(fileID, fmt, 'HeaderLines', 1);
  fclose(fileID);
  % Interpret the time strings and convert it to int64 (tt2000).
  offTime = EpochTT(cell2mat(C{1}));
  % Offset start time (based on ROI).
  time1 = offTime.ttns;
  % Let each offset be valid until 5 us before next offset begin
  % NOTE: While it could be as small as 1ns with int64 representation it
  % can in reallity not be that small as "interp1" in Matlab R2013b req.
  % double representation.
  % (Highest MMS burst rate used as of 2017/05/05 is 16384Hz => dt=61 us).
  time2 = time1(2:end) - int64(5000);
  % Let the last offset be valid "forever" after (interp1 with "linear"
  % and "extrap" require one extra datapoint to ensure it is interpolated
  % as a static value and not a linear trend between the penultimate and
  % ultimate offset value.) Add one extra point one year after the last.
  time2(end+1) = time2(end) + int64(365*86400e9);
  data3 = double([C{2}==1; C{2}==1]); % Repeated data, Look for when it is equal to "1", indicating Solar Wind.
  [timeSort, indSort] = sort([time1; time2]); % Almost repeated time (5 us diff), then sorted
  [sw.time, indUniq] = unique(timeSort); % Ensure no duplicated values
  dataSort = data3(indSort, :); % Sorted data (based on time)
  sw.wake = dataSort(indUniq, :); % Ensure no duplicated values (based on time)

  % Check to see if we are reprocessing any old times, before our "regions"
  % calibration file starts. In that case fall back to the old enabled
  % times for wake removal.
  if sw.time(1) > tInp(1)
    % SW.wake = 1; <-- run sw_wake code, from time and "onwards"
    % SW.wake = 0; <-- don't run sw_wake code, from time and "onwards"
    % Note: SDC run on an old Matlab, so add a duplicate time and wake
    % indicator at the end of segment (on or off), 1 second before
    % change...
    sw.time = []; sw.wake=[];
    irf.log('warning', 'Fall back to static times as regions file starts after our interval');
    sw.time(1)     = EpochTT('2015-03-13T00:00:00.000000000Z').ttns; sw.wake(1)     = 1; % Before start of mission
    sw.time(end+1) = EpochTT('2017-04-19T23:59:59.000000000Z').ttns; sw.wake(end+1) = 1;   % (enabled to a second before change)
    sw.time(end+1) = EpochTT('2017-04-20T00:00:00.000000000Z').ttns; sw.wake(end+1) = 0; % Disable, (orbit is such that no S/W interactions are expected)
    sw.time(end+1) = EpochTT('2017-08-31T18:59:59.000000000Z').ttns; sw.wake(end+1) = 0;   % (disabled to a second before change)
    sw.time(end+1) = EpochTT('2017-08-31T19:00:00.000000000Z').ttns; sw.wake(end+1) = 1; % Enable for one hour (S/W interactions seen in data).
    sw.time(end+1) = EpochTT('2017-08-31T19:59:59.000000000Z').ttns; sw.wake(end+1) = 1;   % (enabled to a second before change)
    sw.time(end+1) = EpochTT('2017-08-31T20:00:00.000000000Z').ttns; sw.wake(end+1) = 0; % Disabled until further notice.
    sw.time(end+1) = EpochTT('2057-08-31T00:00:00.000000000Z').ttns; sw.wake(end+1) = sw.wake(end); % (repeat last entry until end of time...)
  end

else
  % Did not locate any file, use old hard coded values
  switch scId
    case {1,2,3,4}
      % SW.wake = 1; <-- run sw_wake code, from time and "onwards"
      % SW.wake = 0; <-- don't run sw_wake code, from time and "onwards"
      % Note: SDC run on an old Matlab, so add a duplicate time and wake
      % indicator at the end of segment (on or off), 1 second before
      % change...
      irf.log('warning', 'Fall back to static times as no regions file was found.');
      sw.time(1)     = EpochTT('2015-03-13T00:00:00.000000000Z').ttns; sw.wake(1)     = 1; % Before start of mission
      sw.time(end+1) = EpochTT('2017-04-19T23:59:59.000000000Z').ttns; sw.wake(end+1) = 1;   % (enabled to a second before change)
      sw.time(end+1) = EpochTT('2017-04-20T00:00:00.000000000Z').ttns; sw.wake(end+1) = 0; % Disable, (orbit is such that no S/W interactions are expected)
      sw.time(end+1) = EpochTT('2017-08-31T18:59:59.000000000Z').ttns; sw.wake(end+1) = 0;   % (disabled to a second before change)
      sw.time(end+1) = EpochTT('2017-08-31T19:00:00.000000000Z').ttns; sw.wake(end+1) = 1; % Enable for one hour (S/W interactions seen in data).
      sw.time(end+1) = EpochTT('2017-08-31T19:59:59.000000000Z').ttns; sw.wake(end+1) = 1;   % (enabled to a second before change)
      sw.time(end+1) = EpochTT('2017-08-31T20:00:00.000000000Z').ttns; sw.wake(end+1) = 0; % Disabled until further notice.

      sw.time(end+1) = EpochTT('2057-08-31T00:00:00.000000000Z').ttns; sw.wake(end+1) = sw.wake(end); % (repeat last entry until end of time...)
    otherwise
      error('Invalid scId')
  end
end

tInp = double(tInp - sw.time(1));
swtime = double(sw.time - sw.time(1));

indSW = interp1(swtime, sw.wake, tInp, 'previous', 'extrap');

indSW = logical(indSW>0);

end
