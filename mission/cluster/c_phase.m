function [res, tInts] = c_phase(time,phase_2)
%C_PHASE  compute spin phase from PHASE_2
%
% [PHA, T_INTS] = C_PHASE(TIME,PHASE_2)
%
% NaNs are returned when phase cannot be computed
% T_INTS contains time intervals with constant spin rate
%
% See also: MMS_DEFATT_PHASE

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

%% Constants
SPIN_PERIOD_MAX = 4.5; % std 4.5. on commissioning year 2000 value 6.0
SPIN_PERIOD_MIN = 3.6; % std 3.6. on commissioning year 2000 value 2.0
SPIN_PERIOD_NOMINAL = 4; % rpm
DT_MAX = 2*SPIN_PERIOD_MAX; % allow extrapolation for max DT_MAX seconds
SPIN_GAP_MAX = 1200; % max gap in phase we tolerate
MAX_SPIN_PERIOD_CHANGE = SPIN_PERIOD_NOMINAL*0.001;
ERR_PHA_MAX = 0.2; % Error in phase (deg) from fitting

%% Prepare
phaInp = [];
verify_input();
tInts = []; % Output
t0 = phaInp(1,1); tPhase_2 = phaInp(:,1) - t0; targetTime = time - t0;
tPhase_2(tPhase_2<targetTime(1)-DT_MAX/2|tPhase_2>targetTime(end)+DT_MAX/2)=[];
flagSpinRateStable = false; angleError = [];

%% Main loop
tStart = targetTime(1); tStep = targetTime(end) - tStart;
phaseOut = zeros(size(targetTime))*NaN; prevIntExtended = false;
spinPeriodLast = []; iLastOkPoint = []; spinPeriod = [];
while tStart<targetTime(end)
  tStop = tStart + tStep;
  iPhaTmp = tPhase_2>=tStart-DT_MAX/2 & tPhase_2<tStop+DT_MAX/2;
  tPhaTmp = tPhase_2(iPhaTmp);
  if length(tPhaTmp)<=1, tStart = tStop; continue; end
  
  if isempty(iLastOkPoint), iOutTmp = targetTime <= tStop;
  else, iOutTmp = targetTime<=tStop & targetTime>targetTime(iLastOkPoint);
  end
  if ~any(iOutTmp)
    tStart = tStop; tStep = targetTime(end) - tStart;
    continue;
  end
  
  gaps = find(diff(tPhaTmp)>SPIN_GAP_MAX, 1);
  if ~isempty(gaps)
    irf_log('proc','Gaps in phase_2, cannot compute phase'),
    res = []; return
  end
  irf_log('proc',sprintf('Processing %d points (%s -- %s)',...
    length(tPhaTmp),epoch2iso(tStart+t0), epoch2iso(tStop+t0)))
  comp_spin_rate()
  if ~flagSpinRateStable || length(tPhaTmp)>2 && ~isempty(spinPeriodLast) &&...
      abs(spinPeriod-spinPeriodLast) > MAX_SPIN_PERIOD_CHANGE
    if tStep > SPIN_PERIOD_MAX, tStep = tStep/2; continue % Reduce the time step
    else, irf_log('proc','Advancing to next data point')
    end
  else % All good
    % Extend the interval with times of stable spin
    if tStop<targetTime(end) && ~prevIntExtended
      sAEOld = std(angleError); mAEOld = median(angleError);
      idxStart = find(tPhase_2>=tStart-DT_MAX/2,1,'first');
      angleError = comp_angle_error(tPhase_2(idxStart:end));
      iSpinup = find(abs(angleError-mAEOld)>3*sAEOld,1,'first');
      if ~isempty(iSpinup) && iSpinup>3
        tStepTmp = tPhase_2(idxStart+iSpinup-3)-tStart;
        if tStepTmp>tStep
          irf_log('proc',...
            sprintf('Extending the interval by %.1f sec',tStepTmp-tStep))
          tStep = tStepTmp; prevIntExtended = true;
          continue
        end
      end
    end
    prevIntExtended = false;
    polyfit_phase()
    irf_log('proc',sprintf('Spin period %.4f sec (%s -- %s)',...
      spinPeriod,epoch2iso(tStart+t0), epoch2iso(tStop+t0)))
    tInts = [tInts; tStart+t0 tStop+t0]; %#ok<AGROW>
  end
  iLastOkPoint = find(iOutTmp,1,'last');
  tStart = tStop; tStep = targetTime(end) - tStart;
end % Main loop
res = [time phaseOut];

%% Help functions
  function comp_spin_rate()
    flagSpinRateStable = false;
    while true
      comp_spin_period()
      if isempty(spinPeriod)
        irf_log('proc','Cannot determine spin period!'), return
      end
      angleError = comp_angle_error(tPhaTmp);
      medianAngErr = median(angleError); stdAngErr = std(angleError);
      if stdAngErr>ERR_PHA_MAX
        irf_log('proc',...
          sprintf('Std(angleError): %.4f \n',stdAngErr))
        return;
      elseif medianAngErr<ERR_PHA_MAX, flagSpinRateStable = true; return,
      end
      irf_log('proc',...
        sprintf('Median(angleError): %.4f \n',medianAngErr))
      iOut = [];
      find_outliers()
      if isempty(iOut), return, end
      tPhaTmp(iOut) = [];
      if length(tPhaTmp)<2, return, end
    end
    function find_outliers()
      iOut = find( abs(angleError-medianAngErr) > 3*stdAngErr );
      if isempty(iOut), return, end
      irf_log('proc',sprintf('Found %d outliers [TOTAL]',length(iOut)))
      iiIn = find(diff(diff(iOut))==0);
      iOut = setdiff(iOut,iOut(unique([iiIn  iiIn+1 iiIn+2])));
      irf_log('proc',...
        sprintf('Found %d outliers [SINGLE and DOUBLE]',length(iOut)))
    end
    function comp_spin_period()
      dd = diff(tPhaTmp); med = median(dd);
      dd(abs(dd-med)>0.01*med) = []; % remove outliers
      spinPeriod = mean(dd);
      if isnan(spinPeriod) || ...
          spinPeriod > SPIN_PERIOD_MAX || spinPeriod < SPIN_PERIOD_MIN
        if ~isnan(spinPeriod)
            irf_log('proc',sprintf('Found bad spin period: %f s',spinPeriod));
        end
        spinPeriod = [];
      end
    end
  end
  function angleError = comp_angle_error(pha)
    angleError = mod(360.0*(pha-pha(1))/spinPeriod,360);
    angleError = abs(angleError);
    angleError = min([angleError';360-angleError']);
  end
  function polyfit_phase()
    phaseOut(iOutTmp) = ...
      mod((targetTime(iOutTmp)-tPhaTmp(1))/spinPeriod,1)*360;
    spinPeriodLast = spinPeriod;
  end
  function verify_input()
    if size(time,1)>1 && size(time,2)>1, error('t must be a vector'), end
    time = time(:);
    if size(phase_2,1)<2, error('not enough points in phase_2'), end
    phaInp = phase_2(phase_2(:,2)==0,:);
    if isempty(phaInp), error('not enough points in phase_2'), end
  end
end