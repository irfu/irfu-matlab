function res = c_phase(time,phase_2)
%C_PHASE  compute spin phase from PHASE_2
%
% PHA = C_PHASE(TIME,PHASE_2)
%
% See also: MMS_DEFATT_PHASE

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

%% Constants
SPIN_PERIOD_MAX = 4.3;
SPIN_PERIOD_MIN = 3.6;
SPIN_PERIOD_NOMINAL = 4; % rpm
DT_MAX = 2*SPIN_PERIOD_MAX; % allow extrapolation for max DT_MAX seconds
SPIN_GAP_MAX = 300; % max gap in phase we tolerate
MAX_SPIN_PERIOD_CHANGE = SPIN_PERIOD_NOMINAL*0.001;
ERR_PHA_MAX = 0.5; % Error in phase (deg) from fitting

%% Prepare
verify_input();

t0 = phaInp(1,1); tPhase_2 = phaInp(:,1) - t0; targetTime = time - t0;
tPhase_2(tPhase_2<targetTime(1)-DT_MAX|tPhase_2>targetTime(end)+DT_MAX)=[];

%% Main loop
tStart = targetTime(1); tStep = targetTime(end) - tStart; 
phaseOut = zeros(size(targetTime))*NaN; 
spinPeriodLast = []; iLastOkPoint = []; spinPeriod = [];
while tStart<targetTime(end)
  tStop = tStart + tStep;
  iPhaTmp = tPhase_2>=tStart-DT_MAX & tPhase_2<tStop+DT_MAX;
  tPhaTmp = tPhase_2(iPhaTmp);
  if length(tPhaTmp)<=1, tStart = tStop; continue; end
  
  if isempty(iLastOkPoint), iOutTmp = targetTime <= tStop;
  else iOutTmp = targetTime<=tStop & targetTime>targetTime(iLastOkPoint);
  end
  if ~any(iOutTmp), 
    tStart = tStop; tStep = targetTime(end) - tStart; 
    continue; 
  end
  
  %XXX TODO: add handling of gaps
  gaps = find(diff(tPhaTmp)>SPIN_GAP_MAX, 1);
  if ~isempty(gaps), error('gaps'), end
  
  comp_spin_rate()
  if ~flagSpinRateStable || ~isempty(spinPeriodLast) &&...
      abs(spinPeriod-spinPeriodLast) > MAX_SPIN_PERIOD_CHANGE
    if tStep > SPIN_PERIOD_MAX, tStep = tStep/2; continue % Reduce the time step
    else
      irf_log('proc','Do not know how to proceed')
      error('Do not know how to proceed')
    end  
  else % All good
    polyfit_phase()
  end
  iLastOkPoint = find(iOutTmp,1,'last'); 
  tStart = tStop; tStep = targetTime(end) - tStart;
end % Main loop
res = [time phaseOut];

%% Help functions
  function comp_spin_rate()
    flagSpinRateStable = 1;
    dd = diff(tPhaTmp); med = median(dd); 
    dd(abs(dd-med)>0.01*med) = []; % remove outliers
    spinPeriod = mean(dd);
    if isnan(spinPeriod) || ...
        spinPeriod > SPIN_PERIOD_MAX || spinPeriod < SPIN_PERIOD_MIN
      errS = 'Cannot determine spin period!';
      irf.log('critical',errS), error(errS)
    end
    diffangle = mod(360.0*tPhaTmp/spinPeriod,360);
    diffangle = abs(diffangle);
    diffangle = min([diffangle';360-diffangle']);
    if median(diffangle)>ERR_PHA_MAX, flagSpinRateStable = 0; 
      fprintf('Median diff: %.4f \n',median(diffangle)), 
    end
  end
  function polyfit_phase()
    phaseOut(iOutTmp) = mod(targetTime(iOutTmp)/spinPeriod,1)*360;
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