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
    timeTmp(2,:) = [2015 05 23 17 59 59];
    dac.max(2) = -130; dac.min(2) = -190;
    og.max(2)  = -3.9;  og.min(2) = -4.1;
    ig.max(2)  = -7.9;  ig.min(2) = -8.1;
    timeTmp(3,:) = [2015 05 23 18 00 00];
    dac.max(3) = -115; dac.min(3) = -125;
    og.max(3)  = -3.9;  og.min(3) = -4.1;
    ig.max(3)  = -7.9;  ig.min(3) = -8.1;
    % Bias changed to -110nA on 2015/07/23T15 but not really bad data
    timeTmp(4,:) = [2015 07 23 14 59 59];
    dac.max(4) = -115; dac.min(4) = -125;
    og.max(4)  = -3.9;  og.min(4) = -4.1;
    ig.max(4)  = -7.9;  ig.min(4) = -8.1;
    timeTmp(5,:) = [2015 07 23 15 00 00];
    dac.max(5) = -105; dac.min(5) = -115;
    og.max(5)  = -3.9;  og.min(5) = -4.1;
    ig.max(5)  = -7.9;  ig.min(5) = -8.1;
    % Bias changed to -90nA on 2015/08/03T16 but not really bad data
    timeTmp(6,:) = [2015 08 03 15 59 59];
    dac.max(6) = -105;  dac.min(6) = -115;
    og.max(6)  = -3.9;  og.min(6) = -4.1;
    ig.max(6)  = -7.9;  ig.min(6) = -8.1;
    timeTmp(7,:) = [2015 08 03 16 00 00];
    dac.max(7) = -85;  dac.min(7) = -95;
    og.max(7)  = -3.9;  og.min(7) = -4.1;
    ig.max(7)  = -7.9;  ig.min(7) = -8.1;
    % Bias changed to -110nA on 2016/01/26T18:55:16 but not really bad data
    timeTmp(8,:) = [2016 01 26 18 55 15];
    dac.max(8) = -85;  dac.min(8) = -95;
    og.max(8)  = -3.9;  og.min(8) = -4.1;
    ig.max(8)  = -7.9;  ig.min(8) = -8.1;
    timeTmp(9,:) = [2016 01 26 18 55 16];
    dac.max(9) = -105; dac.min(9) = -115;
    og.max(9) = -3.9;   og.min(9) = -4.1;
    ig.max(9) = -7.9;   ig.min(9) = -8.1;
    % Bias changed to -90nA on 2016/09/07T16:19:11
    timeTmp(10,:) = [2016 09 07 16 19 10];
    dac.max(10) = -105; dac.min(10) = -115;
    og.max(10) = -3.9;   og.min(10) = -4.1;
    ig.max(10) = -7.9;   ig.min(10) = -8.1;
    timeTmp(11,:) = [2016 09 07 16 19 11];
    dac.max(11) = -85; dac.min(11) = -95;
    og.max(11) = -3.9;  og.min(11) = -4.1;
    ig.max(11) = -7.9;  ig.min(11) = -8.1;
    % Bias changed to -80nA on p3 and -90nA on p1,p2,p4 on
    % 2016/10/18T14:48:10
    timeTmp(12,:) = [2016 10 18 14 48 09];
    dac.max(12) = -85; dac.min(12) = -95;
    og.max(12)  = -3.9; og.min(12) = -4.1;
    ig.max(12)  = -7.9; ig.min(12) = -8.1;
    timeTmp(13,:) = [2016 10 18 14 48 10];
    dac.max(13) = -75; dac.min(13) = -95; %p3 -80, p1,2,4 -90
    og.max(13)  = -3.9; og.min(13) = -4.1;
    ig.max(13)  = -7.9; ig.min(13) = -8.1;
    % Bias changed to -70nA on p3 and -80nA on p1,2 and -90 on p4 on
    % 2017/05/24T19:32:01
    timeTmp(14,:) = [2017 05 24 19 32 00];
    dac.max(14) = -75; dac.min(14) = -95;
    og.max(14)  = -3.9; og.min(14) = -4.1;
    ig.max(14)  = -7.9; ig.min(14) = -8.1;
    timeTmp(15,:) = [2017 05 24 19 32 01];
    dac.max(15) = -10; dac.min(15) = -200; %p3 -70, p1,2 -80, p4 -90
    og.max(15)  = -3.9; og.min(15) = -4.1;
    ig.max(15)  = -7.9; ig.min(15) = -8.1;
    timeTmp(16,:) = [2017 06 24 19 32 01];
    dac.max(16) = -10; dac.min(16) = -200;
    og.max(16)  = -3.9; og.min(16) = -4.1;
    ig.max(16)  = -7.9; ig.min(16) = -8.1;

  case 2
    % Define new limits for good bias setting on MMS 2
    % Rationale for changing:
    % Bias changed to -100nA on 2015/06/25T02 due to bad data, the really
    % bad data started already around 2015/05/20T03
    timeTmp(2,:) = [2015 05 20 02 59 59];
    dac.max(2) = -130; dac.min(2) = -190;
    og.max(2)  = -3.9; og.min(2) = -4.1;
    ig.max(2)  = -7.9; ig.min(2) = -8.1;
    timeTmp(3,:) = [2015 05 20 03 00 00];
    dac.max(3) = -95; dac.min(3) = -105;
    og.max(3)  = -3.9; og.min(3) = -4.1;
    ig.max(3)  = -7.9; ig.min(3) = -8.1;
    % Bias changed to -90nA on 2015/07/23T16:30 but not really bad data
    timeTmp(4,:) = [2015 07 23 16 29 59];
    dac.max(4) = -95; dac.min(4) = -105;
    og.max(4)  = -3.9; og.min(4) = -4.1;
    ig.max(4)  = -7.9; ig.min(4) = -8.1;
    timeTmp(5,:) = [2015 07 23 16 30 00];
    dac.max(5) = -85; dac.min(5) = -95;
    og.max(5)  = -3.9; og.min(5) = -4.1;
    ig.max(5)  = -7.9; ig.min(5) = -8.1;
    % Bias changed to -70nA on 2015/08/03T20 but not really bad data
    timeTmp(6,:) = [2015 08 03 19 59 59];
    dac.max(6) = -85; dac.min(6) = -95;
    og.max(6)  = -3.9; og.min(6) = -4.1;
    ig.max(6)  = -7.9; ig.min(6) = -8.1;
    timeTmp(7,:) = [2015 08 03 20 00 00];
    dac.max(7) = -65; dac.min(7) = -75;
    og.max(7)  = -3.9; og.min(7) = -4.1;
    ig.max(7)  = -7.9; ig.min(7) = -8.1;
    % Bias changed to -90nA on 2016/01/26T20:55:36 but not really bad data
    timeTmp(8,:) = [2016 01 26 20 55 35];
    dac.max(8) = -65; dac.min(8) = -75;
    og.max(8)  = -3.9; og.min(8) = -4.1;
    ig.max(8)  = -7.9; ig.min(8) = -8.1;
    timeTmp(9,:) = [2016 01 26 20 55 36];
    dac.max(9) = -85; dac.min(9) = -95;
    og.max(9)  = -3.9; og.min(9) = -4.1;
    ig.max(9)  = -7.9; ig.min(9) = -8.1;
    % Bias changed to -90nA on p1, -70nA on p2,p3,p4 on 2016/09/07T17:35:29
    timeTmp(10,:) = [2016 09 07 17 35 28];
    dac.max(10) = -85; dac.min(10) = -95;
    og.max(10)  = -3.9; og.min(10) = -4.1;
    ig.max(10)  = -7.9; ig.min(10) = -8.1;
    timeTmp(11,:) = [2016 09 07 17 35 29];
    dac.max(11) = -65; dac.min(11) = -95;
    og.max(11)  = -3.9; og.min(11) = -4.1;
    ig.max(11)  = -7.9; ig.min(11) = -8.1;
    % Bias change -90nA on p1, -60nA on p2,p3,p4 on 2017/05/24T22:48:09
    timeTmp(12,:) = [2017 05 24 22 48 08];
    dac.max(12) = -65; dac.min(12) = -95;
    og.max(12)  = -3.9; og.min(12) = -4.1;
    ig.max(12)  = -7.9; ig.min(12) = -8.1;
    timeTmp(13,:) = [2017 05 24 22 48 09]; % p1 -90, p234 -60
    dac.max(13) = -10; dac.min(13) = -200;
    og.max(13)  = -3.9; og.min(13) = -4.1;
    ig.max(13)  = -7.9; ig.min(13) = -8.1;
    timeTmp(14,:) = [2017 06 24 22 48 09];
    dac.max(14) = -10; dac.min(14) = -200;
    og.max(14)  = -3.9; og.min(14) = -4.1;
    ig.max(14)  = -7.9; ig.min(14) = -8.1;

  case 3
    % Define new limits for good bias setting on MMS 3
    % Rationale for changing:
    % Bias changed to -120nA on 2015/06/26T16 due to bad data, the really
    % bad data started already around 2015/06/17T07
    timeTmp(2,:) = [2015 06 17 06 59 59];
    dac.max(2) = -130; dac.min(2) = -190;
    og.max(2)  = -3.9;  og.min(2) = -4.1;
    ig.max(2)  = -7.9;  ig.min(2) = -8.1;
    timeTmp(3,:) = [2015 06 17 07 00 00];
    dac.max(3) = -115; dac.min(3) = -125;
    og.max(3)  = -3.9;  og.min(3) = -4.1;
    ig.max(3)  = -7.9;  ig.min(3) = -8.1;
    % Bias changed to -110nA on 2015/07/24T15 but not really bad data
    timeTmp(4,:) = [2015 07 24 14 59 59];
    dac.max(4) = -115; dac.min(4) = -125;
    og.max(4)  = -3.9;  og.min(4) = -4.1;
    ig.max(4)  = -7.9;  ig.min(4) = -8.1;
    timeTmp(5,:) = [2015 07 24 15 00 00];
    dac.max(5) = -105; dac.min(5) = -115;
    og.max(5)  = -3.9;  og.min(5) = -4.1;
    ig.max(5)  = -7.9;  ig.min(5) = -8.1;
    % Bias changed to -90nA on 2015/08/04T02 but not really bad data
    timeTmp(6,:) = [2015 08 04 01 59 59];
    dac.max(6) = -105; dac.min(6) = -115;
    og.max(6)  = -3.9;  og.min(6) = -4.1;
    ig.max(6)  = -7.9;  ig.min(6) = -8.1;
    timeTmp(7,:) = [2015 08 04 02 00 00];
    dac.max(7) = -85; dac.min(7) = -95;
    og.max(7)  = -3.9; og.min(7) = -4.1;
    ig.max(7)  = -7.9; ig.min(7) = -8.1;
    % Bias changed to -110nA on 2016/01/26T14:48:17 but not really bad data
    timeTmp(8,:) = [2016 01 26 14 48 16];
    dac.max(8) = -85; dac.min(8) = -95;
    og.max(8)  = -3.9; og.min(8) = -4.1;
    ig.max(8)  = -7.9; ig.min(8) = -8.1;
    timeTmp(9,:) = [2016 01 26 14 48 17];
    dac.max(9) = -105; dac.min(9) = -115;
    og.max(9)  = -3.9; og.min(9) = -4.1;
    ig.max(9)  = -7.9; ig.min(9) = -8.1;
    % Bias changed to -90nA on p1 and -110nA on p2,p3,p4 on
    % 2016/09/07T18:52:24
    timeTmp(10,:) = [2016 09 07 18 52 23];
    dac.max(10) = -105; dac.min(10) = -115;
    og.max(10)  = -3.9; og.min(10) = -4.1;
    ig.max(10)  = -7.9; ig.min(10) = -8.1;
    timeTmp(11,:) = [2016 09 07 18 52 24];
    dac.max(11) = -85; dac.min(11) = -115;
    og.max(11)  = -3.9; og.min(11) = -4.1;
    ig.max(11)  = -7.9; ig.min(11) = -8.1;
    % Bias changed to -80nA on p1, -90nA on p2, -100nA on p3, -110nA on p4
    % on 2017/05/25T00:56:16
    timeTmp(12,:) = [2017 05 25 00 56 15];
    dac.max(12) = -85; dac.min(12) = -115;
    og.max(12)  = -3.9; og.min(12) = -4.1;
    ig.max(12)  = -7.9; ig.min(12) = -8.1;
    timeTmp(13,:) = [2017 05 25 00 56 16]; %p1 -80, p2 -90, p3 -100, p4 -110
    dac.max(13) = -10; dac.min(13) = -200;
    og.max(13)  = -3.9; og.min(13) = -4.1;
    ig.max(13)  = -7.9; ig.min(13) = -8.1;
    timeTmp(14,:) = [2017 06 25 00 56 16];
    dac.max(14) = -10; dac.min(14) = -200;
    og.max(14)  = -3.9; og.min(14) = -4.1;
    ig.max(14)  = -7.9; ig.min(14) = -8.1;

  case 4
    % Define new limits for good bias setting on MMS 4
    % Rationale for changing:
    % Bias changed to -140nA on 2015/08/03T15 but not really bad data
    %% CHECK WITH FAST & SLOW data
    timeTmp(2,:) = [2015 08 03 14 59 59];
    dac.max(2) = -130; dac.min(2) = -190;
    og.max(2)  = -3.9;  og.min(2) = -4.1;
    ig.max(2)  = -7.9;  ig.min(2) = -8.1;
    timeTmp(3,:) = [2015 08 03 15 00 00];
    dac.max(3) = -135; dac.min(3) = -145;
    og.max(3)  = -3.9;  og.min(3) = -4.1;
    ig.max(3)  = -7.9;  ig.min(3) = -8.1;
    % Bias changed to -160nA on 2016/01/26T16:39:44 but not really bad data
    timeTmp(4,:) = [2016 01 26 16 39 43];
    dac.max(4) = -135; dac.min(4) = -145;
    og.max(4)  = -3.9;  og.min(4) = -4.1;
    ig.max(4)  = -7.9;  ig.min(4) = -8.1;
    timeTmp(5,:) = [2016 01 26 16 39 44];
    dac.max(5) = -155; dac.min(5) = -165;
    og.max(5)  = -3.9;  og.min(5) = -4.1;
    ig.max(5)  = -7.9;  ig.min(5) = -8.1;
    % Bias changed to -140nA on 2016/09/07T20:05:00
    timeTmp(6,:) = [2016 09 07 20 04 59];
    dac.max(6) = -155; dac.min(6) = -165;
    og.max(6)  = -3.9;  og.min(6) = -4.1;
    ig.max(6)  = -7.9;  ig.min(6) = -8.1;
    timeTmp(7,:) = [2016 09 07 20 05 00];
    dac.max(7) = -135; dac.min(7) = -145;
    og.max(7)  = -3.9;  og.min(7) = -4.1;
    ig.max(7)  = -7.9;  ig.min(7) = -8.1;
    % Bias changed to -130nA on 2017/05/25T03:00:55
    timeTmp(8,:) = [2017 05 25 03 00 54];
    dac.max(8) = -135; dac.min(8) = -145;
    og.max(8)  = -3.9;  og.min(8) = -4.1;
    ig.max(8)  = -7.9;  ig.min(8) = -8.1;
    timeTmp(9,:) = [2017 05 25 03 00 55]; % p1234 -130
    dac.max(9) = -10; dac.min(9) = -200;
    og.max(9)  = -3.9;  og.min(9) = -4.1;
    ig.max(9)  = -7.9;  ig.min(9) = -8.1;
    timeTmp(10,:) = [2017 06 25 03 00 55];
    dac.max(10) = -10; dac.min(10) = -200;
    og.max(10)  = -3.9;  og.min(10) = -4.1;
    ig.max(10)  = -7.9;  ig.min(10) = -8.1;

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
