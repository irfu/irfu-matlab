function limits = mms_sdp_limit_bias(scId)
% MMS_SDP_LIMIT_BIAS   get MMS SDP nominal bias limits
%
% Each limit should be applied from from the its start time until further
% notice or a new limit is set. (ie. interp1 with 'previous' and 'extrap').
%
% limits = mms_sdp_limit_bias(scId)
%
% Inputs: scId - MMS spacecraft (1, 2, 3 or 4).
% Output: limits - Struct of TSeries objects
%           .dac - TS with DAC limits [max, min].
%           .og  - TS with OG limits [max, min].
%           .ig  - TS with IG limits [max, min].
%
% See also: TSERIES, MMS_SDP_DMGR/SET_PARAM.

% Verify inputs
narginchk(1,1);

% Limits (pre-launch) used when first turned on. Same for all four s/c.
timeTmp(1,:) = [2015 03 13 00 00 00]; % Launch day
% DAC: 27626TM = -130nA and 25252TM = -190nA. Usually running at -160nA (26439TM).
dac.max(1) = 27625;
dac.min(1) = 25252;
% Outer Guard: 32768TM = 0V and 0TM = -10.4V. Usually running at -4V (20170TM).
og.max(1) = 32767;
og.min(1) = 0;
% Inner Guard: 32768TM = 0V and 0TM = -10.4V. Usually running at -8V (7571TM).
ig.max(1) = 32767;
ig.min(1) = 0;

% Table of new bias limits for each spacecraft
switch scId
  case 1
    % Define new limits for good bias setting on MMS 1
    % Rationale for changing:
    % Bias changed to 28021TM on 2015/06/26T20 due to bad data, really bad data
    % started around 2015/05/25T09
    % Bias changed to 28417TM on 2015/07/23T15
    timeTmp(2,:) = [2015 05 25 09 00 00];
       dac.max(2) = 29000; dac.min(2) = 27625;
       og.max(2) = 32767;  og.min(2) = 0;
       ig.max(2) = 32767;  ig.min(2) = 0;

  case 2
    % Define new limits for good bias setting on MMS 2
    % Rationale for changing:
    % Bias changed to 28812TM on 2015/06/25T03 due to bad data, really bad data
    % started around 2015/05/18T18
    % Bias changed to 29208TM on 2015/07/23T16:30
    timeTmp(2,:) = [2015 05 18 18 00 00];
       dac.max(2) = 30000; dac.min(2) = 28500;
       og.max(2) = 32767;  og.min(2) = 0;
       ig.max(2) = 32767;  ig.min(2) = 0;

  case 3
    % Define new limits for good bias setting on MMS 3
    % Rationale for changing:
    % Bias changed to 28021TM on 2015/06/26T16 due to bad data, really bad data
    % started around 2015/06/19T01
    timeTmp(2,:) = [2015 06 19 01 00 00];
       dac.max(2) = 28500; dac.min(2) = 27625;
       og.max(2) = 32767;  og.min(2) = 0;
       ig.max(2) = 32767;  ig.min(2) = 0;

  case 4
    % Define new limits for good bias setting on MMS 4
    % Rationale for changing: ... Photo e-, bla bla bla
%     timeTmp(2,:) = [2015 04 13 00 00 00];
%        dac.max(2) = 27625; dac.min(2) = 25252;
%        og.max(2) = 32767;  og.min(2) = 0;
%        ig.max(2) = 32767;  ig.min(2) = 0;

  otherwise
    errStr = 'Invalid scId, only numerical 1, 2, 3 or 4 allowed.';
    irf.log('critical',errStr); error(errStr);
end

% Convert timeTmp to EpochTT.
timeTT = EpochTT(irf_time(timeTmp,'vector6>ttns'));
verify_time(timeTT); % Quick check on time.
verify_bias(dac); % Quick check on bias.
verify_bias(og);
verify_bias(ig);
% Convert limits to uint16 and flip them to row vectors.
dac.max = uint16(dac.max)'; dac.min = uint16(dac.min)';
og.max = uint16(og.max)';   og.min = uint16(og.min)';
ig.max = uint16(ig.max)';   ig.min = uint16(ig.min)';

% Define TSeries
limits.dac = TSeries(timeTT,[dac.max, dac.min]);
limits.dac.name = 'DAC (probe) bias limits ([max, min]) for good data.';
limits.og = TSeries(timeTT,[og.max, og.min]);
limits.og.name = 'InnerGuard bias limits ([max, min]) for good data.';
limits.ig = TSeries(timeTT,[ig.max, ig.min]);
limits.ig.name = 'OuterGuard bias limits ([max, min]) for good data.';


% Help functions to verify sanity of limits.
  function verify_time(timeTT)
    % Verify time is monotonically increasing and not FillVal. This can
    % possibly detect a typo in the table above..
    if( any(diff(timeTT.epoch)<=0) || any(timeTT.epoch==int64(-9223372036854775808)))
      errStr = 'Invalid time. Possibly a typo in table of bias limits.';
      irf.log('critical',errStr); error(errStr);
    end
  end

  function verify_bias(bias)
    % Verify bias is not out of range limits and that it is an integer value.
    % This can possibly detect a typo in the table above.
    if( any(bias.max>65535) || any(bias.min>65535) || ...
        any(bias.max<0) || any(bias.min<0) || ...
        any(bias.max <= bias.min) || ...
        any(bias.max~=fix(bias.max)) || any(bias.min~=fix(bias.min)) )
      errStr = 'Invalid bias. Possibly a typo in the table of bias limits.';
      irf.log('critical',errStr); error(errStr);
    end
  end

end
