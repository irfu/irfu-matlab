function [phase, flag, pulse, period, period_flag] = mms_sdp_phase_2(sps, epoch)
% MMS_SDP_PHASE_2 takes sun pulse data (from HK 101) and
% return an array of spin phases (degrees) for each input time epoch
% (TT2000).
%
% Inputs:
%   sps    struct of HK_101 input (Corresponds to 'DATAC.hk_101')
%      .time      - Timestamp of HK packet in TT2000
%      .sunpulse  - Latest sunpulse in TT2000
%      .sunssps   - Sunpulse indicator
%      .iifsunper - CIDP sun period (in us)
%   epoch  array(N), of TT2000 times for which spin phase is desired.
%
% Outputs:
%   phase  array(N), of phase for each corrsponding epoch, in deg.
%   flag   array(N), of flag for each corresponding epoch, in which:
%           0, Interpolated between consecutive sun pulses or average spin
%              rate used to fill gap(-s) < 4 spins.
%           1, Average spin rate used to fill gap(-s) > 4 spins.
%           2, Average spin rate used for gap(-s) < 4 spins, leads to
%              cumulative phase error > 90 degrees.
%           3, True number of spins to fill gap > 4 spins is questionable.
%              Estimates from initial spin rate and final spin rate differ
%              by more than 0.5 spins, but both round to the same integer
%              number of spins.
%           4, True number of spins to fill gap > 4 spins is questionable.
%              Spin period extrapolated from beginning of gap gives
%              cumulative phase error < 10 degrees, but number of spins
%              desagrees with estimate using end spin period, or end spin
%              period is unavailable.
%           5, True number of spins to fill gap(-s) > 4 spins questionable.
%              Spin period extrapolated from end of gap gives cumulative
%              phase error < 10 degrees, but number of spns disagrees with
%              estimate using beginning spin period, or beginning spin
%              period is unavailable.
%           6, True number of spins to fill gap > 4 spins is questionable.
%              Number of spins to fill gap is average of estimates using
%              beginning and ending of spin periods.
%       10-16, Extrapolated beyond first/last sunpulse. 10 is added to
%              quality flag of period used for extrapolation.
%       20-26, Extrapolated > 3 spins beyond first/last sunpulse. 20 is
%              added to quality flag of period used for extrapolation.
%
% Additional optional output(-s), which could be used to create a "cleaner"
% sun pulse file.
%   pulse   array(N), of unique pulse times in TT2000.
%   period  array(N-1), of periods or average periods, to next pulse in ns.
%   period_flag array(N-1), of flags corresponding to each spin period.
%           With same meaning as phase flag.
%
% This code is based on the IDL function mms_fg_sunpulse2phase originally
% written by Ken Bromund, NASA, 2014 for the MMS mission.
%
% Example:
%  [dcephase, dcephase_flag] = MMS_SDP_PHASE_2(DATAC.k_101, DATAC.dce.time);
%
% See also MMS_SDP_DMGR.

% This code updated with source from Ken's mms_fg_sunpulse2phase, dated
% 2016/02/26T22:07:22 UTC as recived in e-mail 2016/09/22.

% Verify number of inputs and outputs
narginchk(2,2);  %
nargoutchk(2,5); % At least the outputs 'phase' and 'flag'.

% Ken Bromunds description is mostly kept in the following segment, as is:
%
% To calculate phase, we need a unique, monotonically increasing array of
% sun pulses, with corresponding spin periods.
% However, we expect to have complications in the HK 0x101 data:
% 1) The same sun-pulse time will be repeated in consecutive 0x101 packets.
% and/or
% 2) Missing sun pulses. E.g. if spin period is less than the 0x101 cadence
% of 10 s. After I wrote the code to handle this, I learned that during SDP
% deployment, the cadence of 0x101 packets will be decreased to 1 s.
% However, I am leaving the code as-is just in case of error in operations,
% and because the same code handles any possible small data gaps.
%
% In addition, this code anticipates larger data gaps.
%
% We are given the spin period in IIFSUNPER (assuming sunssps is 0), but
% this has several issues:
% 1) Not always available (e.g. when in shadow)
% 2) Based on a clock that exhibits temperature drift, not correct TAI
%    seconds. 25 microsec temperature drift obseved in MRT9b data.
% 3) Phase calculation must be based on the sun pulse times, which have
%    jitter, so even if IIFSUNPER is the correct time, it could not be used
%    with the sun pulse times without introducing phase discontinuities.
% For spin period, I use the actual delta-t between sun-pulse times
% provided by the CIDP. These times are relative to TAI, but do exhibit
% some jitter (observed up to 15 us, but sigma is about 1 or 2 us).
%
% I will use IIFSUNPER only for the phase before the initial point of the
% file and after large data gaps (where phase continuity is not an issue).
% It is also used to validate spin period calculations (ie to validate
% correct number of spins between pulses).
%
% A slight improvement could be made by using IIFSUNPER to generate a
% 'pseudo sun pulse time' to reduce the data gap size by 1, whenever
% sun pulses are missing in the smaller gaps. For now, assume that the spin
% rate is changing slowly enough that it is acceptable to assume a
% constant, average spin rate accross small gaps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outline of Algorithm
% 1) Get a unique set of sun pulses from a unique set of HK packets.
% 2) Find spin period that takes us exactly from one pulse to the next,
%    dealing with data gaps.
% 2a) Small data gaps, where it is safe to assume the spin period has not
%     changed significantly enough to allow for uncertainty in the number
%     of spins during the gap.
% 2b) Large data gaps, where spin period might have changed enough that
%     there is uncertainty in the number of spins during the gap.
% 3) Calculate phase, assuming 0.0 degrees when sunpulse occurs. (As per
% e-mail discussion of 2014/09/23).
% 3b) Correct for the sunpulse sensor being +76 degrees of from the {DCS-X,
% DCS-Z} plane (measured along + DCS-Z axis). Correction; subtract 76
% degrees from the calculated value, as per e-mail of 2015/03/20. This way
% the derived phase could be used interchangable with the one derived from
% DEFATT files, as the angle around positive DCS-Z axis, measured from the
% {DCS-X, DCS-Z} plane.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sm_width = 9; % number of points for median smooth. NOTE: Ken uses 10, but that is even.
% Step 1a: remove potential overlap in hk data packets, remove duplicate
%          packets, based on packet time tag.
% first remove sun pulse times with fillval, padval and other rediculously
% early times:
keep = find(sps.sunpulse > spdfcomputett2000([2000 01 01 00 00 00 0 0 0]));
if length(keep) <= sm_width
  logStr = 'Insufficient sun pulse data to determine phase.';
  irf.log('critical', logStr);
  error(logStr);
end
sps.time = sps.time(keep);
sps.sunpulse = sps.sunpulse(keep);
sps.iifsunper = sps.iifsunper(keep);
sps.sunssps = sps.sunssps(keep);

% Locate un-physical times and remove their corresponding sunpulses
% All sunpulses should be some time in the past or at least measured at the
% same instance as the sun pulse triggered the sun sensor. However due to
% some onboard jitter include a 60 seconds margin (typically jitter is less
% than one second, but avoid false positives).
% This problem is most common after long period of instrument off (such as
% during long eclipse periods).
keep = find(sps.sunpulse <= (sps.time + int64(60e9)));
if length(sps.time) ~= length(keep)
  logStr = ['Found ', num2str(length(sps.time)-length(keep)), ...
    ' un-physical pulses (sunpulse >> measurement timestamp). Removing them.'];
  irf.log('warning', logStr);
  sps.time = sps.time(keep);
  sps.sunpulse = sps.sunpulse(keep);
  sps.iifsunper = sps.iifsunper(keep);
  sps.sunssps = sps.sunssps(keep);
end

[sps.time, srt] = sort(sps.time);
[sps.time, usrt] = unique(sps.time);
pulse = sps.sunpulse(srt(usrt));
iifsunper = int64(sps.iifsunper(srt(usrt)))*1000; % convert to nanosec
sunssps = sps.sunssps(srt(usrt));
% ; iifsunper is technically valid if sunssps eq 0 and if iifsunper NE 0,
% ; but I do some sanity checking, too.
valid_iifsunper = sunssps == 0 & iifsunper < 50e9 & iifsunper > 2e9;

% ; Step 1b: We now have a uniqe set of hk packets, but the same sun pulse
% ;          time will be reported in consecutive packets
[pulse, uidx] = unique(pulse);
iifsunper = iifsunper(uidx);
sunssps = sunssps(uidx);
valid_iifsunper = valid_iifsunper(uidx);

% ; step 2: Find spin period.  The initial result will be the spin period
% ;         times the actual number of spins between received pulses.
period = double(circshift(pulse, -1) - pulse);
periodleft = double(pulse - circshift(pulse, 1));
if valid_iifsunper(1)
  periodleft(1) = iifsunper(1);
else
  periodleft(1) = median(period(1:sm_width+1));
end
period(end) = median(period(end-sm_width+1:end));
% ; toss out bad (doubled) pulses received from S/C
smperiod_threshold = IDL_medianfilter(period, sm_width);
% ; note: median smooth does not remove outliers at the edge!  hence, the next two lines of code:
smperiod_threshold(1:fix(sm_width/2)) = median(period(1:sm_width));
smperiod_threshold(end-fix(sm_width/2)+1:end) = median(period(end-sm_width+1:end));
smperiod_threshold = 0.99*smperiod_threshold;
% ; smperiod_threshold *= 0.01 ;; **** for testing.  this disables removal of doubled pulses.
% ; identify points with nominal spinperiod before or after the pulse
nominal_period = (period > smperiod_threshold) | (periodleft > smperiod_threshold);
keep = find(nominal_period); nkeep = length(keep);

ndouble = numel(period)-nkeep;
if(ndouble>0)
  irf.log('warning', ['Removing ', num2str(ndouble), ' repeated sunpulses.']);
  badp = find(nominal_period ~= 1);
  for i=1:min(ndouble,10)
    logStr = sprintf('Removed outlier sunpulse at: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.',...
      spdfbreakdowntt2000(pulse(badp(i))));
    irf.log('notice', logStr);
  end
  % ; iifsunper on the sunpulse after the removed sunpulse is NOT valid
  if(badp(end)>=numel(valid_iifsunper)), badp(end) = []; end
  valid_iifsunper(badp+1) = 0;
end

% ; we might have missed doubled sunpulses if the second pulse occured
% ; within a fraction of a second (~.2 sec at 99% threshod)
% ; which one is right? I don't know.  toss the second one.
pulse = pulse(keep);
period = period(keep);
periodleft = double(pulse - circshift(pulse, 1));
if(valid_iifsunper(1))
  periodleft(1) = iifsunper(1);
else
  periodleft(1) = median(period(1:sm_width+1));
end

quick_threshold = IDL_medianfilter(period, sm_width);
% ; median smooth does not remove outliers at the edge!
quick_threshold(1:fix(sm_width/2)) = median(period(1:sm_width));
quick_threshold(end-fix(sm_width/2)+1:end) = median(period(end-sm_width+1:end));
quick_threshold = 0.50 * quick_threshold;
period_flag = periodleft < quick_threshold;
quick_double = find(period_flag); ndouble = length(quick_double);

if(ndouble>0)
  % ; flag periods before and after removed sunpulse
  irf.log('notice', [num2str(ndouble), ' closely doubled sunpulses, tossing second pulse. Flagging period before and after.']);
  for i=1:min(ndouble, 10)
    logStr = sprintf('Removed second of closely doubled sunpulse pair at: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.',...
      spdfbreakdowntt2000(pulse(quick_double(i))));
    irf.log('notice', logStr);
  end
  keep2 = find(period_flag ~= 1);

  period_flag(quick_double - 1) = 8;
  period_flag(quick_double + 1) = 8;

  % ; iifsunper on the sunpulse before and after the removed sunpulse could be off by 1%
  if quick_double(end)<length(valid_iifsunper)
    valid_iifsunper(quick_double + 1 ) = 0;
  else
    valid_iifsunper(quick_double(1:end-1)+1) = 0;
  end
  valid_iifsunper(quick_double - 1) = 0;
  keep = keep(keep2);
  period_flag = period_flag(keep2);
  pulse = pulse(keep2);
end

% ; period[i] is the period from pulse[i] to pulse[i+1]
period = diff(pulse);

iifsunper = iifsunper(keep);
sunssps = sunssps(keep);
valid_iifsunper = valid_iifsunper(keep);

if valid_iifsunper(1)
  % ; Q: what if this was rejected by the 'if ndouble then' block above?
  % ; A: valid_iifsunper would be set to 0
  pulse = [pulse(1) - iifsunper(1); pulse];
  period = [iifsunper(1); period];
  period_flag = [0; period_flag];
  iifsunper = [0; iifsunper];
  valid_iifsunper = [0; valid_iifsunper];
  sunssps = [3; sunssps];
end

% ; pulses before or after gaps are more likely to be wrong. best to remove
% ; these with median smooth of entire data set.

% ; Detect large gaps, and separate data into segments between large gaps
medianperiod = median(period);

msmoothperiod = IDL_medianfilter(period, sm_width);

mingap = 4;
% ; a large gap is mingap time the median period
gap = find(period > medianperiod*mingap);
ngap = length(gap);
% ; nsegments = ngap + 1
segs = [0, gap', numel(period)];
% ; insert a pseudo sun pulse time and period based on IIFSUNPER at the
% ; beginning of each segment (except the first one, because that's been done)
% ; but we won't do this if iifsunper doesn't agree with the period of the
% ; pulses in the segment.
% ; THIS IS USED TO HELP FIND THE PERIOD AT THE END OF A GAP
for i=1:ngap
  segstart = segs(i) + 1;
  if valid_iifsunper(segstart) && segstart < segs(end) && ...
      abs(iifsunper(segstart) - msmoothperiod(segstart)) < 0.05*1e9
    % ; insert into pulse, period, and period_flag arrays
    pseudopulse = pulse(segstart) - iifsunper(segstart);
    pulse = [pulse(1:segstart); pseudopulse; pulse(segstart+1:end)];
    period = [period(1:segstart); iifsunper(segstart+1); period(segstart+1:end)];
    period_flag = [period_flag(1:segstart); 0; period_flag(segstart+1:end)];
    % ; insert 0 into iifsunper and valid_iifsunper arrays, 3 into ssps
    iifsunper = [iifsunper(1:segstart); 0; iifsunper(segstart+1:end)];
    valid_iifsunper = [valid_iifsunper(1:segstart); 0; valid_iifsunper(segstart+1:end)];
    sunssps = [sunssps(1:segstart); 3; sunssps(segstart+1:end)];
    % ; and don't forget to insert a fill point into the msmoothperiod array, too,
    % ; so that it stays the same length as the iifsunper array.
    msmoothperiod = [msmoothperiod(1:segstart); NaN; msmoothperiod(segstart+1:end)];
    segs(i+2:end) = segs(i+2:end) + 1;
  end
end

spinno = 0:numel(pulse);

if(ngap > 0)
  logStr = sprintf(['%i gaps > %d sec. (%d spins) found in sun pulse data.',...
    ' Processing semi-continuous segments.'], ngap, medianperiod*mingap/1e9, ...
    mingap);
  irf.log('warning', logStr);
end

% ; step 2a: Process each segment between large gaps separately.
% ;          Use IIFSUNPER, if available, or median smooth of periods
% ;          between pulses to determine number of spins between pulses.
for i = 1:ngap+1
  nseg = segs(i+1)-segs(i)-1;  % number of periods in this segment
  if (nseg <= 0), continue; end

  segind = 1+segs(i)+1:nseg+segs(i)+1; % Index, period data of segment.
  segperiod = period(segind);
  segvalid_iifsunper = valid_iifsunper(segind+1);
  segiifsunper = iifsunper(segind+1);
  segind_pulse = (1:nseg+1) + segs(i) + 1; % Index, pulse data of segment.
  segpulse = pulse(segind_pulse);

  % Median smooth the spin period data (this will reject periods that are
  % spuriously long, due to missed sun pulses).
  mwin = 7; % This should be an odd number!

  if(mwin <= nseg)
    mwin2 = floor(mwin/2);
    % TN: MEDIAN in IDL with second argument width applies a median
    % filter. In Matlab is this done by medfilt1 which is found in signal
    % processing toolbox.
    mperiod = IDL_medianfilter(segperiod, 7);
    % ; Don't trust the first and last half-window points of the results of
    % ; IDL's median smooth.
    mperiod(1:mwin2) = mperiod(mwin2+1);
    mperiod(end-mwin2+1:end) = mperiod(end-mwin2);
  else
    % ; segment is shorter than window, so use median of available points
    mperiod = double(segperiod);
    mperiod(:) = median(segperiod);
  end

  % ; Use IIFSUNPER if it is valid
  use_iif = find(segvalid_iifsunper);
  if( isempty(use_iif) )
    logStr = ['Warning: no valid iifsunper values in segment ', num2str(i)];
    irf.log('warning', logStr);
  else
    mperiod(use_iif) = double(segiifsunper(use_iif));
  end

  % Number of spins between each pulse is time between pulse divided by
  % the best guess at the actual period.
  nspins1 = double(segperiod)./mperiod;

  seg_period_flag = period_flag(segind);

  % ; nspins1 should be nearly an integer value.  If not, issue a warning
  % ; and set a flag.
  spin_warn = find( abs(nspins1 - round(nspins1)) >= 0.25 );
  n_spin_warn = length(spin_warn);
  if(n_spin_warn > 0)
    irf.log('warning','spin rate data may be changing too quickly to uniquely determine phase at: ');
    for mi = 1:min(9, n_spin_warn)
      logStr = sprintf('%04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.',...
        spdfbreakdowntt2000(segpulse(spin_warn(mi))));
      irf.log('warning', logStr);
    end
    if(n_spin_warn > 10)
      logStr = ['Error output is truncated: 10 out of ', ...
        num2str(n_spin_warn), ' times shown.'];
      irf.log('warning', logStr);
    end
    seg_period_flag(spin_warn) = 2;
  end

  %    ; TODO (maybe) begin ;;;;;;;;;;;;;;;;;;;
  %    ; insert a pseudo sun pulse time and period based on IIFSUNPER to
  %    ; fill in missed sun pulse time.
  %    ; pseudo = where(round(nspins1) gt 1 and valid_iifsunper[1:*], n_pseudo)
  %    ; insert into pulse, period, and period_flag arrays
  %    ; insert 0 into iifsunper and valid_iifsunper arrays, 3 into ssps
  %    ; adjust segs array indices to accunt for inserted points
  %    ; re-calculate period[pseudo]
  %    ; decrement nspins1[pseudo]
  %    ; and furter debug the mess that this incremental improvment would cause...
  %    ; TODO (mabye) end   ;;;;;;;;;;;;;;;;;;;;

  % ; calculate the average period between pulses.
  nspins1 = int64(round(nspins1));
  period(segind) = period(segind) ./ nspins1;
  period_flag(segind) = seg_period_flag;
  addspins = find(nspins1 > 1); nadd=length(addspins);
  for gi = 1:nadd
    % ; increase spinno after the gap by the number of missed sunpulses.
    spinno(segind(addspins(gi)):end) = spinno(segind(addspins(gi)):end) + double(nspins1(addspins(gi)));
  end
end

% ; step 2b: it is not safe to apply the median smooth method across large gaps.
% ;          Instead, we check to see if it is reasonable to assume that
% ;          the spin rate has not changed during the gap.

for i = 1:ngap
  gapidx = segs(i+1);
  % ; 1) find period before/after gap
  % ; 2) determine nspins in gap based on those two periods
  %
  % ; 1a) find period at beginning of gap:  median of last 5 spins or
  % ;     IIFSUNPER of last pulse before gap
  if gapidx -1 ~= segs(i)
    period1 = median(period(max(gapidx-1-5,segs(i)+1) : gapidx-1));
    % ; if pulse time at beginning of gap is out of family, then fix it
    if abs(period(gapidx-1)-period1) > 0.05*1e9
      pulse(gapidx) = pulse(gapidx-1) + period1;
      period(gapidx-1) = period1;
      % ;TODO : flag it
    end
  elseif(valid_iifsunper(gapidx))
    period1 = double(iifsunper(gapidx));
  else
    logStr = sprintf(['Processing gap %i, cannot determine spin period at ',...
      'beginning of gap: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'],...
      i, spdfbreakdowntt2000(pulse(gapidx)));
    irf.log('warning', logStr);
    period1 = NaN;
  end

  % ; 1b) find period at end of gap:  median of IIFSUNPER (from pseudopulse)
  % ;     of first pulse after gap, plus next 4 spin periods.
  if gapidx + 1 ~= segs(i+2) && gapidx < length(period)
    period2 = median(period(gapidx+1 : min(gapidx+1+5, segs(i+2)-1) ));
    % ; if pulse time at end of gap is out of family, then fix it
    if abs(period(gapidx+1) - period2) > 0.05*1e9
      pulse(gapidx+1) = pulse(gapidx+2) - period2;
      period(gapidx+1) = period2;
      % ;TODO : flag it
    end
  else
    % ; note: unlike for period1, we don't have an option to take period2
    % ;   from iifsunper, because if it looked valid, we would already have
    % ;   inserted this.
    logStr = sprintf(['Processing gap %i, cannot determine spin period at ',...
      'end of gap: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'],...
      i, spdfbreakdowntt2000(pulse(gapidx+1)));
    irf.log('warning', logStr);
    period2 = NaN;
  end

  period(gapidx) = pulse(gapidx+1) - pulse(gapidx);

  % ; 2) determine number of spins in gap, working forwards and backwards.
  nspins1 = period(gapidx)/period1;
  nspins2 = period(gapidx)/period2;

  if isfinite(nspins1) && isfinite(nspins2) && round(nspins1) == round(nspins2)
    nspins = (nspins1+nspins2)/2;
    period_flag(gapidx) = 1;
    if abs(nspins1 - nspins2) > 0.5
      logStr = sprintf(['Processing gap %i spin rate interpolation ', ...
        'questionable: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'], ...
        i, spdfbreakdowntt2000(pulse(gapidx)));
      irf.log('warning', logStr);
      period_flag(gapidx) = 3;
    end
  elseif isfinite(nspins1) && abs(nspins1 - round(nspins1)) < 0.25*0.111
    % ; if the period from either end of the gap seems correct to 10 degrees,
    % ; then accept it.
    nspins = nspins1;
    logStr = sprintf(['Processing gap %i accepting spin period from ', ...
      'beginning of gap at: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'], ...
      i, spdfbreakdowntt2000(pulse(gapidx)));
    irf.log('warning', logStr);
    period_flag(gapidx) = 4;
  elseif finite(nspins2) && abs(nspins2 - round(nspins2)) < 0.25*0.111
    nspins = nspins2;
    logStr=sprintf(['Processing gap %i accepting spin period from ', ...
      'end of gap at: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'], ...
      i, spdfbreakdowntt2000(pulse(gapidx)));
    irf.log('warning', logStr);
    period_flag(gapidx) = 5;
  else
    nspins = (nspins1+nspins2)/2;
    logStr=sprintf(['Processing gap %i period at end of gap is ', ...
      'incompatible with period at beginning. Using average spin period ', ...
      'as best guess at: %04i-%02i-%02iT%02i:%02i:%02i:%03i.%03i.%03iZ.'], ...
      i, spdfbreakdowntt2000(pulse(gapidx)));
    irf.log('warning', logStr);
    period_flag(gapidx) = 6;
  end

  period(gapidx) = period(gapidx) ./ round(nspins); % average spin period over the interval where there are no sunpulses
  % ; for Hermite interpolation:
  % ; ; increase spinno after the gap by the number of missed sunpulses.
  % ; spinno[gapidx+1:*] += round(nspins-1)
  % ; ; set period at beginning and end of gap to best guess: will be used by hermite interpolation.
  % ; period[gapidx] = period1
  % ; period[gapidx+1] = period2
end

%  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
%  ;step 3:  calculate spin phase
%  ;         This based on the spin period and the time since last sun
%  ;         pulse for each input epoch
%  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

% Transform to 'double'. It is possible to do this earlier but perhaps best
% to keep int64 until all TT2000 specific calls are done (spdfbreakdowntt2000).
epoch = double(epoch);
pulse = double(pulse);
period = double(period);

%per = value_locate(pulse, epoch)
per = interp1(pulse, 1:length(pulse), epoch, 'nearest');

% Pre allocate output
phase = NaN(size(epoch));
flag = int16(epoch);

%  ; phase for BCS is sun pulse phase -76 degrees
% 2015/03/20 Mail from Ken, changing standpoint from 2014/09/23, sunpulse
% timestamp is when sensor see sun. This is not in the {BCS-X, BCS-Z} plane
% but at the +76 degrees of from BCS-X where the sensor is located.
% Therefor, in order to give correct phase, shift the calculated value by
% negative 76 deg. Then remap to interval [0-360) degrees.
dss2bcs = -76;

%  inrange = where(per ge 0 and per lt n_elements(pulse)-1, n_inrange)
inrange = find( ~isnan(per) & per < length(pulse) );
n_inrange = length(inrange);

if isempty(inrange)
  % Empty overlap meaning HK 101 sunpulses do not cover any of the
  % EDP data, may even be hours off... This was the case for
  % mms4_edp_fast_l1b_dce_20170904_v1.4.0.cdf (covering 00.00-00.30 UTC)
  % while the corresponding HK 101 file
  % mms4_fields_hk_l1b_101_20170904_v0.5.0.cdf
  % covering 08:25:21 to 08:44:42 UTC.
  % Interpolation  is not an option for these cases, return error and wait
  % for more instrument data along with more hk data (or final DEFATT).
  logStr = 'Missmatch between hk101 sunpulses and L1b data, no overlap in time.';
  irf.log('critical', logStr);
  error(logStr);
end

if n_inrange > 0

  % ; linear interpolate between each pulse
  phase(inrange) = (torow(epoch(inrange)) - torow(pulse(per(inrange))))./torow(period(per(inrange)))*360 + dss2bcs;
  % ;    ; hermite interpolate between each pulse
  % ;    t0 = pulse[inrange[0]]
  % ;    phase[inrange] = hermite(double(pulse[inrange]-t0), spinno[inrange]*360.d, double(epoch[inrange]-t0), $
  % ;      fderiv=360./period[inrange]) + dss2bcs
  % ;    store_data, 'phase_h', data={x:time_double(epoch[inrange], /tt2000), y:phase[inrange]}
  flag(inrange) = period_flag(per(inrange));
  phase(inrange) = mod(phase(inrange), 360);
end

% ; epochs before and after first/last sunpulse are special cases
if inrange(1)>1
  % ; the first pulse will be a 'pseudo pulse', with period from IIFSUNPER,
  % ; if available.
  % ; Otherwise, it will be the first pulse from the CDIP, with spin period
  % ; calculated between first and second pulses
  beforefirst = 1:inrange(1)-1;
  period0 = period(1);
  period_flag0 = 10 + period_flag(1); % Add initial 10 to flag.
  phase(beforefirst) = (epoch(beforefirst)-pulse(1))/period0*360 + dss2bcs;
  flag(beforefirst) = period_flag0;
  warn = find(phase < -360*3);
  if( ~isempty(warn) )
    logStr='Warning: extrapolating more than 3 spin periods before first sun pulse.';
    irf.log('warning', logStr);
    flag(warn) = flag(warn) + 10; % Add another 10 to flag (atleast 10 has already been added).
  end
  phase(beforefirst) = phase(beforefirst)+ceil(-min(phase/360))*360;
  % phase in the interval 0-360 deg.
  phase(beforefirst) = mod(phase(beforefirst), 360);
end

if( inrange(end) < length(epoch) )
  afterlast = inrange(end)+1:length(epoch);
  % Some missing times at the end
  % To calculate phase after last pulse, use spin period calculated
  % between last and second to last pulses (which will actually be
  % IIFSUNPER if there was a gap before the last pulse)
  period0 = period(end);
  flag(afterlast) = 10 + period_flag(end);
  phase(afterlast) = (epoch(afterlast)-pulse(end))/period0*360 + dss2bcs;
  warn = find(phase > 360*3);
  if( ~isempty(warn) )
    logStr='Warning: extrapolating more than 3 spin periods after last sun pulse.';
    irf.log('warning', logStr);
    flag(warn) = flag(warn) + 10;
  end
  % phase in the interval 0-360 deg.
  phase(afterlast) = mod(phase(afterlast), 360);
end

end

function out = IDL_medianfilter(sigIn, width)
% Local function doing similar thing as IDL median in regards to its
% filtering (argument width) and without mean for even numbers.
% To be run in Matlab. Author: T.Nilsson, IRFU, 2014

% Verify arguments
narginchk(2,2)
if(length(sigIn)<width || width<1 || mod(width,2)==0 || floor(width)~=width)
  % Incorrect input arguments.
  err_str = 'Wrong input(-s) to IDL_MEDIAN';
  irf.log('critical', err_str);
  error('MATLAB:MMS_SDP_PHASE_2:IDL_MEDIAN:INPUT', err_str);
end

% Preallocate output
out = zeros(length(sigIn),1);

width_half = (width-1) / 2;

% Edges, beginning and end, left untouched.
out(1:width_half) = sigIn(1:width_half);
out(end-width_half:end) = sigIn(end-width_half:end);

% Middle bit, median of the sorted values of each segment of size width.
for ii = 1+width_half:length(sigIn)-width_half
  localSeg = sort(sigIn(ii-width_half:ii + width_half));
  out(ii) = median(localSeg);
end

end
