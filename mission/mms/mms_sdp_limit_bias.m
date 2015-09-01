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
% DAC: Usually running at -160nA (26439TM).
dac.max(1) = -130;  dac.min(1) = -190;  % nA
% Outer Guard: Usually running at -4V (20170TM).
og.max(1) = -3.9;    og.min(1) = -4.1;  % V
% Inner Guard: Usually running at -8V (7571TM).
ig.max(1) = -7.9;    ig.min(1) = -8.1;  % V

% Table of new bias limits for each spacecraft
switch scId
  case 1
    % Define new limits for good bias setting on MMS 1
    % Rationale for changing:
    % Bias changed to -120nA on 2015/06/26T20 due to bad data, the really
    % bad data started already around 2015/05/23T18
    timeTmp(2,:) = [2015 05 23 18 00 00];
       dac.max(2) = -115; dac.min(2) = -125;
       og.max(2)  = -3.9;  og.min(2) = -4.1;
       ig.max(2)  = -7.9;  ig.min(2) = -8.1;
    % Bias changed to -110nA on 2015/07/23T15 but not really bad data
    timeTmp(3,:) = [2015 07 23 15 00 00];
       dac.max(3) = -105; dac.min(3) = -115;
       og.max(3)  = -3.9;  og.min(3) = -4.1;
       ig.max(3)  = -7.9;  ig.min(3) = -8.1;
    % Bias changed to -90nA on 2015/08/03T16 but not really bad data
    timeTmp(4,:) = [2015 08 03 16 00 00];
       dac.max(4) = -85;  dac.min(4) = -95;
       og.max(4)  = -3.9;  og.min(4) = -4.1;
       ig.max(4)  = -7.9;  ig.min(4) = -8.1;

  case 2
    % Define new limits for good bias setting on MMS 2
    % Rationale for changing:
    % Bias changed to -100nA on 2015/06/25T02 due to bad data, the really
    % bad data started already around 2015/05/20T03
    timeTmp(2,:) = [2015 05 20 03 00 00];
       dac.max(2) = -95; dac.min(2) = -105;
       og.max(2)  = -3.9; og.min(2) = -4.1;
       ig.max(2)  = -7.9; ig.min(2) = -8.1;
    % Bias changed to -90nA on 2015/07/23T16:30 but not really bad data
    timeTmp(3,:) = [2015 07 23 16 30 00];
       dac.max(3) = -85; dac.min(3) = -95;
       og.max(3)  = -3.9; og.min(3) = -4.1;
       ig.max(3)  = -7.9; ig.min(3) = -8.1;
    % Bias changed to -70nA on 2015/08/03T20 but not really bad data
    timeTmp(4,:) = [2015 08 03 20 00 00];
       dac.max(4) = -65; dac.min(4) = -75;
       og.max(4)  = -3.9; og.min(4) = -4.1;
       ig.max(4)  = -7.9; ig.min(4) = -8.1;

  case 3
    % Define new limits for good bias setting on MMS 3
    % Rationale for changing:
    % Bias changed to -120nA on 2015/06/26T16 due to bad data, the really
    % bad data started already around 2015/06/17T07
    timeTmp(2,:) = [2015 06 17 07 00 00];
       dac.max(2) = -115; dac.min(2) = -125;
       og.max(2)  = -3.9;  og.min(2) = -4.1;
       ig.max(2)  = -7.9;  ig.min(2) = -8.1;
    % Bias changed to -110nA on 2015/07/24T15 but not really bad data
    timeTmp(3,:) = [2015 07 24 15 00 00];
       dac.max(3) = -105; dac.min(3) = -115;
       og.max(3)  = -3.9;  og.min(3) = -4.1;
       ig.max(3)  = -7.9;  ig.min(3) = -8.1;
    % Bias changed to -90nA on 2015/08/04T02 but not really bad data
    timeTmp(4,:) = [2015 08 04 02 00 00];
       dac.max(4) = -85; dac.min(4) = -95;
       og.max(4)  = -3.9; og.min(4) = -4.1;
       ig.max(4)  = -7.9; ig.min(4) = -8.1;

  case 4
    % Define new limits for good bias setting on MMS 4
    % Rationale for changing:
    % Bias changed to -140nA on 2015/08/03T15 but not really bad data
    %% CHECK WITH FAST & SLOW data
    timeTmp(2,:) = [2015 08 03 15 00 00];
        dac.max(2) = -135; dac.min(2) = -145;
        og.max(2)  = -3.9;  og.min(2) = -4.1;
        ig.max(2)  = -7.9;  ig.min(2) = -8.1;

  otherwise
    errStr = 'Invalid scId, only numerical 1, 2, 3 or 4 allowed.';
    irf.log('critical',errStr); error(errStr);
end

% Convert timeTmp to EpochTT.
timeTT = EpochTT(irf_time(timeTmp,'vector6>ttns'));
verify_time(timeTT); % Quick check on time.
verify_bias_I(dac); % Quick check on bias.
verify_bias_V(og);
verify_bias_V(ig);
% Convert limits to single and flip them to row vectors.
dac.max = single(dac.max)'; dac.min = single(dac.min)';
og.max  = single(og.max)';   og.min = single(og.min)';
ig.max  = single(ig.max)';   ig.min = single(ig.min)';

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

  function verify_bias_V(bias)
    % Verify bias is not out of range limits. This can possibly detect a
    % typo in the table above.
    if( any(bias.max>10.403) || any(bias.min>10.403) || ...
        any(bias.max<-10.403) || any(bias.min<-10.403) || ...
        any(bias.max <= bias.min) )
      errStr = 'Invalid bias. Possibly a typo in the table of bias limits.';
      irf.log('critical',errStr); error(errStr);
    end
  end

  function verify_bias_I(bias)
    % Verify bias is not out of range limits. This can possibly detect a
    % typo in the table above. Limits are from the usable range, as defined
    % in "ADP/SDP BEB Command Interface v1.7"
   if( any(bias.max>110) || any(bias.min>110) || ...
        any(bias.max<-550) || any(bias.min<-550) || ...
        any(bias.max <= bias.min) )
      errStr = 'Invalid bias. Possibly a typo in the table of bias limits.';
      irf.log('critical',errStr); error(errStr);
    end
  end

end
