function out = irf_spin_epoch(data,phase, fCut)
%IRF_SPIN_EPOCH  compute a spin epoch
%
%  Out = irf_spin_epoch(Data,Phase, fCut)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DATA    - tseries of data (vector, scalar)
%          PHASE   - phase corresponding to DATA time. 
%          FCUT    - frequency for high-pass filter

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<3, fCut = []; end

STEPS_PER_DEG=10;
N_SPINS_MEDIAN = 31; % number of spins for moving median

phaUnw = unwrap(double(phase.data)*pi/180)*180/pi;
fxPha = ((phaUnw(1)-rem(phaUnw(1),360)):1/STEPS_PER_DEG:...
      (phaUnw(end)-rem(phaUnw(end),360)+360))';
fxPha(end) = [];
epoch0 = data.time(1); epochTmp = double(data.time-epoch0);
tFxPha = interp1(phaUnw,epochTmp,fxPha);
nSpins = length(fxPha)/360/STEPS_PER_DEG;

out = data; 
if ~isempty(fCut), out = out.filt(fCut,0,[],5); end
eRes = interp1(epochTmp,double(out.data),tFxPha);
idxOK = ~isnan(tFxPha); nComps = size(eRes,2);
modelOut = zeros(length(epochTmp),nComps);
for iComp = 1:nComps
  eResM = reshape(eRes(:,iComp),360*STEPS_PER_DEG,nSpins);
  model = movmedian(eResM,N_SPINS_MEDIAN,2,'omitnan');  %FIXME: movmedian() is not part of Matlab R2013a installed at SDC!!
  model = reshape(model,numel(model),1);
  modelOut(:,iComp) = interp1(tFxPha(idxOK),model(idxOK),epochTmp,'linear','extrap');
end
out.data = modelOut;
end