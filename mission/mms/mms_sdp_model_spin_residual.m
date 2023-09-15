function [ModelOut,Model360] = mms_sdp_model_spin_residual(Dce,Dcv,Phase,signals,sampleRate)
%MMS_SDP_MODEL_SPIN_RESIDUAL  create a spin residual
%
%  modelOut = mms_sdp_model_spin_residual(dce,dcv,phase,signals,sampleRate)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DCE     - structure with fields time, e12, e34
%          DCV     - structure with fields time, v1, v2, v3, v4
%          PHASE   - phase corresponding to DCE time.
%          SIGNALS - cell array with list of signals to proceess,
%                    e.g {'e12', 'e34'} or {'v1', 'v3'}

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

MMS_CONST = mms_constants;

%% compute model
epoch0 = Dce.time(1); epochTmp = double(Dce.time-epoch0);

% Spinfits setup
MAX_IT = 3;      % Maximum of iterations to run fit
N_TERMS = 3;     % Number of terms to fit, Y = A + B*sin(wt) + C*cos(wt) +..., must be odd.
MIN_FRAC = 0.20; % Minumum fraction of points required for one fit (minPts = minFraction * fitInterv [s] * samplerate [smpl/s] )
FIT_EVERY = 5*10^9;   % Fit every X nanoseconds.
FIT_INTERV = double(MMS_CONST.Limit.SPINFIT_INTERV); % Fit over X nanoseconds interval.

sdpPr = {'v1', 'v2','v3','v4'};
Sfit = struct(sdpPr{1}, [],sdpPr{2}, [], sdpPr{3}, [],sdpPr{4}, []);

% Calculate minumum number of points req. for one fit covering fitInterv
minPts = MIN_FRAC * sampleRate * FIT_INTERV/10^9; % "/10^9" as fitInterv is in [ns].

% Calculate first timestamp of spinfits to be after start of dce time
% and evenly divisable with fitEvery.
% I.e. if fitEvery = 5 s, then spinfit timestamps would be like
% [00.00.00; 00.00.05; 00.00.10; 00.00.15;] etc.
% For this one must rely on spdfbreakdowntt2000 as the TT2000 (int64)
% includes things like leap seconds.
t1 = spdfbreakdowntt2000(Dcv.time(1)); % Start time in format [YYYY MM DD HH MM ss mm uu nn]
% Evenly divisable timestamp with fitEvery after t1, in ns.
t2 = ceil((t1(6)*10^9+t1(7)*10^6+t1(8)*10^3+t1(9))/FIT_EVERY)*FIT_EVERY;
% Note; spdfcomputett2000 can handle any column greater than expected,
% ie "62 seconds" are re-calculated to "1 minute and 2 sec".
t3.sec = floor(t2/10^9);
t3.ms  = floor((t2-t3.sec*10^9)/10^6);
t3.us  = floor((t2-t3.sec*10^9-t3.ms*10^6)/10^3);
t3.ns  = floor(t2-t3.sec*10^9-t3.ms*10^6-t3.us*10^3);
% Compute what TT2000 time that corresponds to, using spdfcomputeTT2000.
t0 = spdfcomputett2000([t1(1) t1(2) t1(3) t1(4) t1(5) t3.sec t3.ms t3.us t3.ns]);

phaRad = unwrap(Phase.data*pi/180);

STEPS_PER_DEG = 1; phaShift=STEPS_PER_DEG/2;
phaDegUnw = phaRad*180/pi;
phaFixed = (fix(phaDegUnw(1)):STEPS_PER_DEG:fix(phaDegUnw(end)))' ...
  + phaShift;
pha360 = 0:STEPS_PER_DEG:360; pha360 = pha360' + phaShift; pha360(end) = [];
n360 = length(pha360);
timeTmp = interp1(phaDegUnw,epochTmp,phaFixed);
phaFixed(isnan(timeTmp)) = []; timeTmp(isnan(timeTmp)) = [];
phaFixedWrp = mod(phaFixed,360);


for signal = signals
  sig = signal{:};
  if isempty(intersect(signal,{'e12','e34','v1','v2','v3','v4'}))
    errS = ['invalid signal: ' sig]; irf.log('critical',errS), error(errS)
  end

  phaseRadTmp = phaRad;
  if( (Dcv.time(1)<=t0) && (t0<=Dcv.time(end)))
    bits = MMS_CONST.Bitmask.SWEEP_DATA;
    if sig(1)=='e'
      dataIn = Dce.(sig).data;
      dataIn = mask_bits(dataIn, Dce.(sig).bitmask, bits);
      timeIn = Dce.time;
    else
      dataIn = Dcv.(sig).data;
      dataIn = mask_bits(dataIn, Dcv.(sig).bitmask, bits);
      timeIn = Dcv.time;
    end
    idxGood = ~isnan(dataIn);
    [tSfit, Sfit.(sig), ~, ~, ~] = ...
      mms_spinfit_m(MAX_IT, minPts, N_TERMS, double(timeIn(idxGood)),...
      double(dataIn(idxGood)), phaseRadTmp(idxGood), FIT_EVERY, FIT_INTERV, t0);
  else
    warnStr = sprintf(['Too short time series:'...
      ' no data cover first spinfit timestamp (t0=%i)'],t0);
    irf.log('warning', warnStr);
  end

  sfitR = interp1(double(tSfit-epoch0),Sfit.(sig),epochTmp);
  spinFitComponent = sfitR(:,1) + sfitR(:,2).*cos(phaRad) + sfitR(:,3).*sin(phaRad);

  Model360.(sig) = zeros(n360,1);
  spinRes = interp1(epochTmp,double(double(dataIn))-spinFitComponent,timeTmp);
  for idx=1:n360
    dataTmp = spinRes(phaFixedWrp==pha360(idx));
    Model360.(sig)(idx) = median(dataTmp(~isnan(dataTmp)));
  end
  %Interpolate through the ADP spike for V1..4
  if sig(1)=='v'
    angProbeShadow = mod(MMS_CONST.Phaseshift.(['p' sig(2)])*180/pi+180,360);
    idxAdp = find(pha360<angProbeShadow+1 & pha360>angProbeShadow-1); % +/- 1 deg
    idxOK = [idxAdp(1)-1 idxAdp(end)+1];
    Model360.(sig)(idxAdp) = interp1(pha360(idxOK),Model360.(sig)(idxOK),...
      pha360(idxAdp));
  end
  ModelOut.(sig) = interp1([-phaShift; pha360; 360+phaShift],...
    [Model360.(sig)(end); Model360.(sig); Model360.(sig)(1)]',Phase.data);
end

end