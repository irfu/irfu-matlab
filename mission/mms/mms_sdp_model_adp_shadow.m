function modelOut = mms_sdp_model_adp_shadow(dce,phase,signals)
%MMS_SDP_MODEL_ADP_SHADOW  create a model for ADP shadow
%
%  modelOut = mms_sdp_model_adp_shadow(dce,phase,signals)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DCE     - structure with fields time, e12, e34,
%                    or time, v1, v2, v3 and v4.
%          PHASE   - phase corresponding to DCE time.
%          SIGNALS - cell array with list of signals to proceess,
%                    e.g {'e12', 'e34'} or {'e12', 'p123'} or
%                    {'v1','v2','v3','v4'}

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

STEPS_PER_DEG=10;
N_SPINS_MEDIAN = 31; % number of spins for moving median
DETR_DPHA = 2; % window for detrending (+/- from the center) in deg
SHA_DPHA = .9; % window for shadow model (+/- from the center) in deg
SHA_OFF = 0.1; % offset shadow position (empirical) in deg

phaUnw = unwrap(phase.data*pi/180)*180/pi;
fxPha = ((phaUnw(1)-rem(phaUnw(1),360)):1/STEPS_PER_DEG:...
  (phaUnw(end)-rem(phaUnw(end),360)+360))';
fxPha(end) = [];
epoch0 = dce.time(1); epochTmp = double(dce.time-epoch0);
tFxPha = interp1(phaUnw,epochTmp,fxPha);
nSpins = length(fxPha)/360/STEPS_PER_DEG;

for iSig = 1:length(signals)
  sig = signals{iSig};
  expShadow = getExpShadow();
  eRes = interp1(epochTmp,double(dce.(sig).data),tFxPha);
  eResM = reshape(eRes,360*STEPS_PER_DEG,nSpins);
  model = zeros(size(eResM));
  phaOne = (0:1/STEPS_PER_DEG:360)'; phaOne(end) = [];
  for sha = expShadow
    % Detrend the data
    idxDetr = (phaOne>sha-DETR_DPHA & phaOne<sha-SHA_DPHA) | ...
      (phaOne<sha+DETR_DPHA & phaOne>sha+SHA_DPHA);
    A = phaOne(idxDetr); X = [ones(length(A),1) A]; x = X\eResM(idxDetr,:);
    idxDetr1 = (phaOne>sha-DETR_DPHA & phaOne<sha+DETR_DPHA);
    eResM(idxDetr1,:) = eResM(idxDetr1,:) - ...
      repmat(x(1,:),length(find(idxDetr1>0)),1) - phaOne(idxDetr1)*x(2,:);
    % Compute moving median
    iModel = (phaOne>sha-SHA_DPHA & phaOne<sha+SHA_DPHA);
    model(iModel,:) = movmedian(eResM(iModel,:),N_SPINS_MEDIAN,2,'omitnan');
  end
  model = reshape(model,numel(model),1);
  idxOK = ~isnan(tFxPha);
  modelOut.(sig) = interp1(tFxPha(idxOK),model(idxOK),epochTmp,'linear','extrap');
end

  function expShadow = getExpShadow()
    % return expected shadow in degrees
    global MMS_CONST
    if isempty(MMS_CONST), MMS_CONST = mms_constants; end
    switch sig
      % Phaseshift.pX = 0 rad when probe X is sunward, i.e. in shade pi rad later
      case 'e12', expShadow = [MMS_CONST.Phaseshift.p1, MMS_CONST.Phaseshift.p2] + pi;
      case 'e34', expShadow = [MMS_CONST.Phaseshift.p3, MMS_CONST.Phaseshift.p4] + pi;
      case 'p123'
        expShadow = [MMS_CONST.Phaseshift.p1, MMS_CONST.Phaseshift.p2, MMS_CONST.Phaseshift.p3] + pi; % p4 lost
        sig = 'e34';
      case 'p124'
        expShadow = [MMS_CONST.Phaseshift.p1, MMS_CONST.Phaseshift.p2, MMS_CONST.Phaseshift.p4] + pi; % p3 lost
        sig = 'e34';
      case 'p134'
        expShadow = [MMS_CONST.Phaseshift.p1, MMS_CONST.Phaseshift.p3, MMS_CONST.Phaseshift.p4] + pi; % p2 lost
        sig = 'e12';
      case 'p234'
        expShadow = [MMS_CONST.Phaseshift.p2, MMS_CONST.Phaseshift.p3, MMS_CONST.Phaseshift.p4] + pi; % p1 lost
        sig = 'e12';
      case 'v1', expShadow = MMS_CONST.Phaseshift.p1 + pi; % p1
      case 'v2', expShadow = MMS_CONST.Phaseshift.p2 + pi; % p2
      case 'v3', expShadow = MMS_CONST.Phaseshift.p3 + pi; % p3
      case 'v4', expShadow = MMS_CONST.Phaseshift.p4 + pi; % p4
      otherwise
        errS = 'unrecognized SIG';
        irf.log('critical',errS), error(errS)
    end
    expShadow = mod(180*expShadow/pi, 360) +SHA_OFF;
  end %getExpShadow()
end
