function res = mms_defatt_phase(defatt,time)
%MMS_DEFATT_PHASE  compute spin phase from DEFATT
%
% PHA = MMS_DEFATT_PHASE(DEFATT)

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
if time(1) > EpochTT2000('2015-06-18T00:00:00.000000Z').epoch
  SPIN_RATE_MAX = MMS_CONST.Spinrate.max;
elseif time(1) > EpochTT2000('2015-01-01T00:00:00.000000Z').epoch
  SPIN_RATE_MAX = MMS_CONST.Spinrate.max_comm; 
else SPIN_RATE_MAX = MMS_CONST.Spinrate.max_deploy;
end
SPIN_RATE_NOMINAL = 3.1; % rpm
DT_MAX = 2*60/SPIN_RATE_NOMINAL; % allow extrapolation for max DT_MAX seconds
STEP_NSPINS_DEF = 30;
MAX_SPIN_RATE_CHANGE = SPIN_RATE_NOMINAL*0.001;

verify_input();

t0 = defatt.time(1); targetTime = double(time-t0)*1e-9;
phaseOut = zeros(size(targetTime))*NaN;
tDefatt = double(defatt.time-t0)*1e-9; phaseDefatt = defatt.zphase;

iOut = ( tDefatt<targetTime(1)-DT_MAX | tDefatt>targetTime(end)+DT_MAX );
tDefatt(iOut) = []; phaseDefatt(iOut) = []; 

spinRateLast = []; iLastGoodPoint = [];
tStep = STEP_NSPINS_DEF*60/SPIN_RATE_NOMINAL; tStart = 0;
while tStart<=targetTime(end)
  tStop = tStart+tStep/2;
  iPhaTmp = tDefatt>=tStart-tStep/2 & tDefatt<tStop;
  tPhaTmp = tDefatt(iPhaTmp); phaTmp = phaseDefatt(iPhaTmp);
  
  if length(tPhaTmp)<=1, tStart = tStop; continue; end
  
  gaps = find(diff(tPhaTmp)>60/SPIN_RATE_MAX);
  if ~isempty(gaps), error('gaps'), end
  
  comp_spin_rate();
  if ~flagSpinRateStable
    error('spinup')
  end
  if ~isempty(spinRateLast) &&...
      abs(spinRate-spinRateLast) > MAX_SPIN_RATE_CHANGE
    error('slow spinup')
  end
  
  % Polyfit
  if isempty(iLastGoodPoint), iOutTmp = targetTime < tStop;
  else iOutTmp = targetTime<tStop & targetTime>targetTime(iLastGoodPoint);
  end
  phaseOut(iOutTmp) = polyval(fitCoef, targetTime(iOutTmp));
  phaseOut(iOutTmp) = mod(phaseOut(iOutTmp),360);
  iLastGoodPoint = find(iOutTmp,1,'last'); tStart = tStop;
end

res = TSeries(EpochTT2000(time),phaseOut);

  function comp_spin_rate()
    flagSpinRateStable = 1;
    
    phc = unwrap(phaTmp*pi/180)*180/pi; fitCoef = polyfit(tPhaTmp,phc,1);
    if isnan(fitCoef(1))
      irf_log('proc','Cannot determine spin period!'), return
    end
    diffangle = mod(phaTmp - polyval(fitCoef,tPhaTmp),360);
    diffangle = abs(diffangle);
    diffangle = min([diffangle';360-diffangle']);
    if max(diffangle)>0.001, flagSpinRateStable = 0; end
    spinRate  = 60/fitCoef(1);
  end
  function verify_input()
    if ~isa(time,'int64'), 
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