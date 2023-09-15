function res = mms_defatt_phase(defatt,time)
%MMS_DEFATT_PHASE  compute spin phase from DEFATT
%
% PHA = MMS_DEFATT_PHASE(DEFATT,TIME)
%
% Returns TSeries

%% Constants
global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
if time(1) > EpochTT('2015-06-18T00:00:00.000000Z').epoch
  SPIN_RATE_MAX = MMS_CONST.Spinrate.max;
elseif time(1) > EpochTT('2015-01-01T00:00:00.000000Z').epoch
  SPIN_RATE_MAX = MMS_CONST.Spinrate.max_comm;
else, SPIN_RATE_MAX = MMS_CONST.Spinrate.max_deploy;
end
SPIN_RATE_NOMINAL = 3.1; % rpm
DT_MAX = 60/SPIN_RATE_NOMINAL; % allow extrapolation for max DT_MAX seconds
STEP_NSPINS_DEF = 30;  TSTEP_MAX = STEP_NSPINS_DEF*60/SPIN_RATE_NOMINAL;
MAX_SPIN_RATE_CHANGE = SPIN_RATE_NOMINAL*0.001;
ERR_PHA_MAX = 0.05; % Error in phase (deg) from fitting
flagSpinRateStable = 1; spinRate = SPIN_RATE_NOMINAL;
%% Prepare
verify_input();

t0 = defatt.time(1);
targetTime = double(time-t0)*1e-9; tStart = targetTime(1);
tDefatt = double(defatt.time-t0)*1e-9; phaseDefatt = defatt.zphase;

phaseOut = zeros(size(targetTime))*NaN;
iOut = ( tDefatt<targetTime(1)-DT_MAX | tDefatt>targetTime(end)+DT_MAX );
tDefatt(iOut) = []; phaseDefatt(iOut) = [];

%% Main loop
spinRateLast = []; iLastOkPoint = []; fitCoef = []; tStep = TSTEP_MAX;
while tStart<=targetTime(end)
  tStop = tStart + tStep;
  iPhaTmp = tDefatt>=tStart-tStep/2 & tDefatt<tStart+tStep*3/2;
  tPhaTmp = tDefatt(iPhaTmp); phaTmp = phaseDefatt(iPhaTmp);
  if length(tPhaTmp)<=1, tStart = tStop; continue; end
  phaTmpUnwrapped = unwrap(phaTmp*pi/180)*180/pi;

  if isempty(iLastOkPoint), iOutTmp = targetTime < tStop;
  else, iOutTmp = targetTime<tStop & targetTime>targetTime(iLastOkPoint);
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
      interp_phase()
    end
  else % All good
    polyfit_phase()
  end
  iLastOkPoint = find(iOutTmp,1,'last'); tStep = TSTEP_MAX; tStart = tStop;
end % Main loop
res = TSeries(EpochTT(time),phaseOut);

%% Help functions
  function comp_spin_rate()
    flagSpinRateStable = 1;
    fitCoef = polyfit(tPhaTmp,phaTmpUnwrapped,1);
    if isnan(fitCoef(1))
      errS = 'Cannot determine spin period!';
      irf.log('critical',errS), error(errS)
    end
    diffangle = mod(phaTmp - polyval(fitCoef,tPhaTmp),360);
    diffangle = abs(diffangle);
    diffangle = min([diffangle';360-diffangle']);
    if median(diffangle)>ERR_PHA_MAX, flagSpinRateStable = 0;
      %fprintf('Median diff: %.4f \n',median(diffangle)),
    end
    spinRate  = 60/fitCoef(1);
  end
  function interp_phase()
    %disp('    >>>>>>>    interpolating >>>>>   -----')
    iOutTmp = targetTime<tStop & targetTime>=tStart;
    phaseOut(iOutTmp) = interp1(tPhaTmp,phaTmpUnwrapped,...
      targetTime(iOutTmp),'linear','extrap');
    phaseOut(iOutTmp) = mod(phaseOut(iOutTmp),360);
    spinRateLast = [];
  end
  function polyfit_phase()
    phaseOut(iOutTmp) = polyval(fitCoef, targetTime(iOutTmp));
    phaseOut(iOutTmp) = mod(phaseOut(iOutTmp),360);
    spinRateLast = spinRate;
  end
  function verify_input()
    if ~isa(time,'int64')
      errStr = 'TIME must be int64 (tt2000)';
      irf.log('critical',errStr), error(errStr)
    end
    if ~isstruct(defatt) || ~isfield (defatt,'time') || ...
        ~isfield (defatt,'zphase') || ~isa(defatt.time,'int64') || ...
        length(defatt.time) ~= length(defatt.zphase)
      errStr = 'DEFATT is not in expected format';
      irf.log('critical',errStr), error(errStr)
    end
    if time(1)>defatt.time(end) || time(end)<defatt.time(1)
      errStr = 'DEFATT and TIME do not overlap';
      irf.log('critical',errStr), error(errStr)
    end
  end
end