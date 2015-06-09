function res = c_phase_new(time,phase_2)
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
SPIN_RATE_MAX = 4.3;
%SPIN_RATE_MIN = 3.6;
SPIN_RATE_NOMINAL = 15; % rpm
DT_MAX = SPIN_RATE_MAX; % allow extrapolation for max DT_MAX seconds
STEP_NSPINS_DEF = 300;  TSTEP_MAX = STEP_NSPINS_DEF*60/SPIN_RATE_NOMINAL;
MAX_SPIN_RATE_CHANGE = SPIN_RATE_NOMINAL*0.001;
ERR_PHA_MAX = 0.5; % Error in phase (deg) from fitting

%% Prepare
verify_input();

phaInp = phase_2(phase_2(:,2)==0,:);
if isempty(phaInp), error('not enough points in phase_2'), end

t0 = phaInp(1,1); 
targetTime = time-t0; tStart = targetTime(1);
tPhase_2 = phaInp(:,1)-t0;

phaseOut = zeros(size(targetTime))*NaN;
iOut = ( tPhase_2<targetTime(1)-DT_MAX | tPhase_2>targetTime(end)+DT_MAX );
tPhase_2(iOut) = [];

%% Main loop
spinRateLast = []; iLastOkPoint = []; spinPeriod = []; tStep = TSTEP_MAX; 
while tStart<=targetTime(end)
  tStop = tStart + tStep;
  iPhaTmp = tPhase_2>=tStart-tStep/2 & tPhase_2<tStart+tStep*3/2;
  tPhaTmp = tPhase_2(iPhaTmp);
  if length(tPhaTmp)<=1, tStart = tStop; continue; end
  
  if isempty(iLastOkPoint), iOutTmp = targetTime < tStop;
  else iOutTmp = targetTime<tStop & targetTime>targetTime(iLastOkPoint);
  end
  if ~any(iOutTmp), tStart = tStop; continue; end
  
  %XXX TODO: add handling of gaps
  gaps = find(diff(tPhaTmp)>60/SPIN_RATE_MAX, 1);
  if ~isempty(gaps), error('gaps'), end
  
  comp_spin_rate()
  if ~flagSpinRateStable || ~isempty(spinRateLast) &&...
      abs(spinRate-spinRateLast) > MAX_SPIN_RATE_CHANGE
    if tStep > DT_MAX, tStep = tStep/2; continue % Reduce the time step
    else
      irf_log('proc','Do not know how to proceed')
      error('Do not know how to proceed')
    end  
  else % All good
    polyfit_phase()
  end
  iLastOkPoint = find(iOutTmp,1,'last'); tStep = TSTEP_MAX; tStart = tStop;
end % Main loop
res = [time phaseOut];

%% Help functions
  function comp_spin_rate()
    flagSpinRateStable = 1;
    dd = diff(tPhaTmp); dd(abs(dd-median(dd))>3*std(dd)) = [];
    spinPeriod = mean(dd);
    if isnan(spinPeriod(1))
      errS = 'Cannot determine spin period!';
      irf.log('critical',errS), error(errS)
    end
    diffangle = mod(360.0*tPhaTmp/spinPeriod,360);
    diffangle = abs(diffangle);
    diffangle = min([diffangle';360-diffangle']);
    if median(diffangle)>ERR_PHA_MAX, flagSpinRateStable = 0; 
      fprintf('Median diff: %.4f \n',median(diffangle)), 
    end
    spinRate  = 60/spinPeriod(1);
  end
  function polyfit_phase()
    phaseOut(iOutTmp) = mod(targetTime(iOutTmp)/spinPeriod,1)*360;
    spinRateLast = spinRate;
  end
  function verify_input()
    if size(time,1)>1 && size(time,2)>1, error('t must be a vector'), end
    if size(phase_2,1)<2, error('not enough points in phase_2'), end
  end
end