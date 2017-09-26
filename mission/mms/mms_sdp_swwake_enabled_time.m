function indSW = mms_sdp_swwake_enabled_time(tInp, scId)
% Returns index "indSW" of when mms S/W wake removal code should run, and
% when it should not.

narginchk(2,2)
if isa(tInp,'GenericTimeArray'), tInp = EpochTT(tInp).ttns; end 

% Table of boom lengths
switch scId
  case {1,2,3,4}
    % SW.wake = 1; <-- run sw_wake code, from time and "onwards"
    % SW.wake = 0; <-- don't run sw_wake code, from time and "onwards"
    % Note: SDC run on an old Matlab, so add a duplicate time and wake
    % indicator at the end of segment (on or off), 1 second before
    % change...
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

tInp = double(tInp - sw.time(1));
swtime = double(sw.time - sw.time(1));

indSW = interp1(swtime, sw.wake, tInp, 'previous', 'extrap');

indSW = logical(indSW>0);

end