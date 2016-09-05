function modelOut = mms_sdp_model_adp_shadow(dce,phase,signals)
%MMS_SDP_MODEL_ADP_SHADOW  create a model for ADP shadow
%
%  modelOut = mms_sdp_model_adp_shadow(dce,phase,signals)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DCE     - structure with fields time, e12, e34
%          PHASE   - phase corresponding to DCE time. 
%          SIGNALS - cell array with list of signals to proceess, 
%                    e.g {'e12', 'e34'} or {'e12', 'p123'}

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
global MMS_CONST
STEPS_PER_DEG=10;
STEP_SPINS = 30;

phaseUnw = unwrap(phase.data*pi/180)*180/pi;
phaseTmp = (fix(phaseUnw(1)):1/STEPS_PER_DEG:ceil(phaseUnw(end)))';
epoch0 = dce.time(1); epochTmp = double(dce.time-epoch0);
timeTmp = interp1(phaseUnw,epochTmp,phaseTmp);
phaseTmp(isnan(timeTmp)) = []; timeTmp(isnan(timeTmp)) = []; 
phaseTmpWrp = mod(phaseTmp,360);

for iSig = 1:length(signals)
  sig = signals{iSig};
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
    otherwise,
      errS = 'unrecognized SIG';
      irf.log('critical',errS), error(errS)
  end
  expShadow = mod(180*expShadow/pi, 360);
  eRes = interp1(epochTmp,double(dce.(sig).data),timeTmp);
  model = zeros(size(phaseTmp));
  
  nData = length(phaseTmpWrp); 
  winPoints = STEP_SPINS*3600*STEPS_PER_DEG;
  stepPoints = fix(winPoints/3); % we advance by 1/3 of the window
  for iS = 1:length(expShadow)
    sha = expShadow(iS);
    idxFirst = 1;
    while true
      idxLast = idxFirst + winPoints -1;
      if idxLast>nData
        idxLast = nData; idxFirst = idxLast - winPoints;
        if idxFirst<1, idxFirst = 1; end
      end
      idxWin = (1:nData)'; idxWin = idxWin>=idxFirst & idxWin<=idxLast;
      idxAssign = (1:nData)';
      idxAssFirst = idxFirst +stepPoints; idxAssLast = idxLast -stepPoints;
      if idxFirst==1, idxAssFirst = 1; end
      if idxLast == nData, idxAssLast = nData; end
      idxAssign = idxAssign>=idxAssFirst & idxAssLast<=idxLast;
      
      idx = idxWin & ( (phaseTmpWrp>sha-2 & phaseTmpWrp<sha-1) | ...
        (phaseTmpWrp>sha+1 & phaseTmpWrp<sha+2) );
      off =  median(eRes(idx));
      for iPha=(-10:1:10)
        pha = sha+iPha/STEPS_PER_DEG;
        idxPha = int16(phaseTmpWrp*10)==int16(pha*10);
        wakeVal = median(eRes(idxWin & idxPha))-off;
        model(idxAssign & idxPha) = wakeVal;
        %if iPha==0
        %  sprintf('Sig %s, sha=%d spin #%.1f amp = %.2f mV/m',sig,sha,...
        %    idxFirst/winPoints,wakeVal)
        %end
      end
      idxFirst = idxFirst+stepPoints;
      if idxFirst>=nData-winPoints+stepPoints, break, end
    end
  end
  modelOut.(sig) = interp1(timeTmp,model,epochTmp,'spline');
end