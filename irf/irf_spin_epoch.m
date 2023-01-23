function out = irf_spin_epoch(varargin)
%IRF_SPIN_EPOCH  compute a spin epoch
%
%  Out = irf_spin_epoch(Data,Phase,'fCut',fCut,'samplefreq',samplefreq)
%
%  Created a model to a disturbace signal caused by the ADP shadow by
%  looking at many spins.
%
%  Input : DATA    - tseries of data (vector, scalar)
%          PHASE   - phase corresponding to DATA time.
%
%  Options:
%          FCUT    - frequency for high-pass filter
%          NSPINS  - number of spins used to construct spin epoch (default is 31)
%          SAMPLEFREQ - sample frequency

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

N_SPINS_MEDIAN = 31; % number of spins for moving median
fCut = []; % Default, no high-pass filter
SAMPLEFREQ = []; % Default

if (nargin < 2)
  nargin
  help irf_spin_epoch;
  return;
end

data=varargin{1};
phase=varargin{2};
args=varargin(3:end);

if numel(args)>0
  flag_have_options=1;
else
  flag_have_options=0;
end

while flag_have_options
  l = 2;
  switch(lower(args{1}))
    case 'fcut'
      if numel(args)>1 && isnumeric(args{2})
        fCut = args{2};
        irf.log('notice',['fCut is set to ' num2str(fCut)]);
      end
    case 'nspins'
      if numel(args)>1 && isnumeric(args{2})
        N_SPINS_MEDIAN = args{2};
        irf.log('notice',['Number of spins is set to ' num2str(N_SPINS_MEDIAN)]);
      end
    case 'samplefreq'
      if numel(args)>1 && isnumeric(args{2})
        SAMPLEFREQ = args{2};
        irf.log('notice',['samplefreq is set to ' num2str(SAMPLEFREQ)]);
      end
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}]);
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args)
    flag_have_options=0;
  end
end


STEPS_PER_DEG=2;

phaUnw = unwrap(double(phase.data)*pi/180)*180/pi;
% Gaps may not have been unwrapped correctly so add a complete revolution
% (360 deg) to ensure unique phase values
idxGap = find(diff(phaUnw)<0);
for iGap=1:length(idxGap)
  phaUnw(idxGap(iGap)+1:end) = phaUnw(idxGap(iGap)+1:end) + 360;
end
fxPha = ((phaUnw(1)-rem(phaUnw(1),360)):1/STEPS_PER_DEG:...
  (phaUnw(end)-rem(phaUnw(end),360)+360))';
fxPha(end) = [];
epoch0 = data.time(1); epochTmp = double(data.time-epoch0);
tFxPha = interp1(phaUnw,epochTmp,fxPha);
nSpins = length(fxPha)/360/STEPS_PER_DEG;

out = data;
if ~isempty(fCut), out = out.filt(fCut,0,SAMPLEFREQ,5); end
eRes = interp1(epochTmp,double(out.data),tFxPha);
idxOK = ~isnan(tFxPha); nComps = size(eRes,2);
modelOut = zeros(length(epochTmp),nComps);
for iComp = 1:nComps
  eResM = reshape(eRes(:,iComp),360*STEPS_PER_DEG,nSpins);
  model = movmedian(eResM,N_SPINS_MEDIAN,2,'omitnan');
  model = reshape(model,numel(model),1);
  modelOut(:,iComp) = interp1(tFxPha(idxOK),model(idxOK),epochTmp,'linear','extrap');
end
out.data = modelOut;
end
