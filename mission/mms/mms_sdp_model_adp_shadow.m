function modelOut = mms_sdp_model_adp_shadow(dce,phase)
%MMS_SDP_MODEL_ADP_SHADOW  create a model for ADP shadow
%
%  modelOut = mms_sdp_model_adp_shadow(dce,phase)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DCE - structure with fields time, e12, e34
%        PHASE - phase corresponding to DCE time. 

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

STEPS_PER_DEG=10;
STEP_SPINS = 30;

phaseUnw = unwrap(phase.data*pi/180)*180/pi;
phaseTmp = (fix(phaseUnw(1)):1/STEPS_PER_DEG:ceil(phaseUnw(end)))';
epoch0 = dce.time(1); epochTmp = double(dce.time-epoch0);
timeTmp = interp1(phaseUnw,epochTmp,phaseTmp);
phaseTmp(isnan(timeTmp)) = []; timeTmp(isnan(timeTmp)) = []; 
phaseTmpWrp = mod(phaseTmp,360);

signals = {'e12','e34'}; expShadow = 150 + [0 180];
for iSig = 1:length(signals)
  if iSig==2, expShadow = expShadow - 90; end % p34
  sig = signals{iSig};
  eRes = interp1(epochTmp,double(dce.(sig).data),timeTmp);
  model = zeros(size(phaseTmp));
  
  nData = length(phaseTmpWrp); 
  winPoints = STEP_SPINS*3600*STEPS_PER_DEG;
  stepPoints = fix(winPoints/3); % we advance by 1/3 of the window
  for iS = 1:length(expShadow)
    sha = expShadow(iS);
    idxFirst = 1;
    while idxFirst<nData-winPoints+stepPoints
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
    end
  end
  modelOut.(sig) = interp1(timeTmp,model,epochTmp,'spline');
end