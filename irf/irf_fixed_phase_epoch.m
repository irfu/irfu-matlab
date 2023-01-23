function res = irf_fixed_phase_epoch(Phase,pha0,phaStep,phaOffs)
%IRF_FIXED_PHASE_EPOCH  create timeseries of constant phase
%
% PHASEFIXED = IRF_FIXED_PHASE_EPOCH(PHASE[,PHA_0,PHA_STEP,PHA_OFFS])
%
% Create a correponding timeline of fixed phase, e.g. 0 to 360 deg,
% which can be used for phase epoch analysis.
%
% Input:
%    PHASE - phase in degrees
%    PHA_0 - reference phase (default 0)
%    PHA_STEP - step in degrees (degault 1)
%    PHA_OFFS - offset from PHA_0 (default PHA_STEP/2)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


if nargin<2, pha0 = 0; end
if nargin<3, phaStep = 1; end
if nargin<4, phaOffs=phaStep/2;
else
  if phaOffs>phaStep
    errS = 'offset>step'; irf.log('critical',errS), error(errS);
  end
end
if isa(Phase,'TSeries'), flagTS = true; else, flagTS = false; end

if flagTS
  phaData = double(Phase.data);
  epoch0 = Phase.time(1);
  epochTmp = Phase.time.epochUnix - epoch0.epochUnix;
else
  epoch0 = Phase(1,1); epochTmp = Phase(:,1) - epoch0;
  if isnan(epoch0)
    errS = 'epoch0 is NaN'; irf.log('critical',errS), error(errS);
  end
  phaData = Phase(:,2);
end

phaDegUnw = unwrap(phaData*pi/180)*180/pi;

% Handle long gaps
N_GAP = 10; % Number of missing points to bec considered a data gap
idxGap = find(diff(epochTmp)>N_GAP*median(diff(epochTmp)));
idxStart = [1; idxGap+1];
idxEnd = [idxGap; length(epochTmp)];

phaFixedAll = []; timeAll = [];
for idx = 1:length(idxStart)
  irf.log('notice',sprintf('Phase gap #%d',idx-1));
  phaFixed =...
    ((round((phaDegUnw(idxStart(idx))-(pha0+phaOffs))/phaStep)*phaStep):...
    phaStep:...
    (round((phaDegUnw(idxEnd(idx))-(pha0+phaOffs))/phaStep)*phaStep))'...
    +pha0+phaOffs;
  idxTmp = idxStart(idx):idxEnd(idx);
  timeTmp = interp1(phaDegUnw(idxTmp),epochTmp(idxTmp),phaFixed);
  phaFixed(isnan(timeTmp)) = []; timeTmp(isnan(timeTmp)) = [];
  timeAll = [timeAll; timeTmp]; phaFixedAll = [phaFixedAll; phaFixed]; %#ok<AGROW>
end
phaFixedWrp = mod(phaFixedAll,360);
timePhaFixed = epoch0 + timeAll;

if flagTS, res = irf.ts_scalar(timePhaFixed,phaFixedWrp);
else, res = [timePhaFixed, phaFixedWrp];
end

end