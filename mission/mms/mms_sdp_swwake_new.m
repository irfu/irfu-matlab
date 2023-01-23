function [wakeModelOut, n_corrected, wakedesc] = mms_sdp_swwake_new(e, pair, phaseDeg, timeIn, NPOINTS,swFlag)
% MMS_SDP_SWWAKE_NEW Correct SDP data for wake in the solar wind
%
% [data, n_corrected, wakedesc] = mms_sdp_swwake_new(e, pair, phase_2, timeIn, sampleRate, solarWindFlag)
%
% Input:
%   e          - raw E-field data
%   pair       - probe pair ('e12' or 'e34')
%   phase_2    - spinphase corresponding to measurements
%   timeIn     - time of measurement (int64, tt2000 ns)
%   sampleRate - nominal sample rate of data (Hz)
%   solarWindFlag - flag(1/0) indicating solar wind per datapoint (e)
%
% Output:
%   data - corrected data
%   n_corrected - numer of spins corrected for wake
%   wakedesc - wake description,
%              as matrix of columns: [time, phase, amplitude, half width]
%
% Wakes are identified by max derivative. First we find a narrow proxy
% wake, correct for it and find a proxy DC field (ground tone). Then we
% subtract the ground tone from the original data and find a final fit for
% the wake. The procedure is performed on five spins, with the resulting
% fit being applied to the spin in the middle.
%
% This program was written in order to improve quality of the sdp data in
% the solar wind.
%
% Based on code for Cluster EFW but re-written for MMS SDP.
% See also: c_efw_swwake

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Original idea by Anders Eriksson.
% Many useful suggestions by Per-Arne Lindqvist.
% Re-written for MMS by Thomas Nilsson.

% ThoNi NOTES:
% Resampling data to resolution 1 deg, sorting data from multiple spins in
% bins of same phase. For MMS however it should be possible to use 0.5
% degrees resolution instead of 1 deg, as MMS Fast mode data have 32 Hz and
% one spin takes about 19-20 second (approx phase resolution = 360/(19x32)
% = 0.59 deg per sample).

narginchk(4,6)
global MMS_CONST;
if(isempty(MMS_CONST)), MMS_CONST = mms_constants; end

wakeModelOut = [];
n_corrected = 0;
wakedesc = [];

% N_EMPTY 0.75 means that we use only spins with more then 75% points.
%N_EMPTY = 0.9; % (Cluster was 0.9)
%MAX_SPIN_PERIOD = 4.3; % sec for Cluster
%MAX_SPIN_PERIOD = 10^9*60/MMS_CONST.Spinrate.min; % = 20 sec
WAKE_MAX_HALFWIDTH = 25; % degrees, (Cluster was 45 deg)
WAKE_MIN_HALFWIDTH = 9; %11;  % degrees, (Cluster was 11 deg)
WAKE_MIN_AMPLITUDE = 0.35; % mV/m, (Cluster was 0.4 mV/m)
WAKE_MAX_AMPLITUDE = 7; % mV/m, (Cluster was 7 mV/m)
if nargin < 5
  NPOINTS = 360; % Number of points per spin
end
if nargin < 6
  swFlag = []; % Number of points per spin
end
NWSPINS = 5; % number of spins we work on

plot_step = 1;
plot_i = 0;
plotflag = false; plotflag_now = false;

switch pair
  case {'e12', 'e34'}
  otherwise
    errStr = 'Pair must be one of: "e12" or "e34"';
    irf.log('critical', errStr);
    error(errStr);
end

% Find numer of spins using the fact that Z-phase is monotonically
% increasing, except for it being modulo 360.
i0 = find(diff(phaseDeg)<0);
if length(i0)<=5
  irf.log('notice', 'Not enough spins for correction.');
  return
end

[dataFixedPha,fixedPha,epochFixedPha,PHASE_OFF] = ...
  mms_interp_fixed_pha(e, timeIn, phaseDeg,360/NPOINTS, pair);

expPhase = 180/pi*MMS_CONST.Phaseshift.(pair) + PHASE_OFF; % XXX PHASE_OFF here is the shift introduce by mms_interp_fixed_pha()
expPhaseIdx = deg2idx(expPhase + WAKE_MAX_HALFWIDTH/2*[-1 1]);
expPhaseIdx = expPhaseIdx(1):expPhaseIdx(2);

if ~isempty(swFlag)
  epoch0 = timeIn(1); epochFixedPhaTmp = double(epochFixedPha-epoch0);
  epochTmp = double(timeIn-epoch0);
  swFlag = interp1(epochTmp,double(swFlag),epochFixedPhaTmp,'linear');
  swFlag = logical(swFlag>0);
end

wakeModel = zeros(size(dataFixedPha));

idx0 = (find(fixedPha==0 | fixedPha==180))';
wakedesc = NaN(length(idx0), 4);

for idx = idx0
  if isnan(dataFixedPha(idx)), continue, end
  if length(wakeModel)-idx<NPOINTS*(NWSPINS)-1, break, end
  ts = epochFixedPha(idx);
  % wakedesc first column is the time stamp of wake.
  wakedesc(find(idx==idx0,1,'first'), 1) = ts;
  tStUTC = irf_time(ts,'ttns>utc');
  
  pha5spins = reshape(fixedPha(idx:(idx+NPOINTS*NWSPINS-1)),NPOINTS,NWSPINS);
  %epoch5spins = reshape(epochFixedPha(idx:(idx+NPOINTS*NWSPINS-1)),NPOINTS,NWSPINS);
  d5spins = reshape(dataFixedPha(idx:(idx+NPOINTS*NWSPINS-1)),NPOINTS,NWSPINS);
  %xxx = fft(d5spins); xxx(5:end-1) = 0;
  %d5spinsAC = d5spins -ifft(xxx, 'symmetric');
  
  idxSpin = (idx:idx+NPOINTS-1) + NPOINTS*2;
  idx1=1:(NPOINTS/2); idx2=(NPOINTS/2+1):NPOINTS;
  idxSpin1 = idxSpin(1) -1 + idx1;
  idxSpin2 = idxSpin(1) -1 + idx2;
  flagSWSpin = isSW(idxSpin);
  
  % Spin in the middle has maximum weigth
  Ki = [.1, .25, .3, .25, .1];
  Ni = Ki;
  if 1 % add more weighting
    av12Prel = sum(d5spins .* repmat(Ki, NPOINTS, 1), 2);
    %Wi = sum(abs(d5spins-repmat(av12Prel,1,5))); % linear weight
    Wi = sum((d5spins-repmat(av12Prel,1,5)).^2); % squared
    Ni([1 2 4 5]) = Ki([1 2 4 5])./Wi([1 2 4 5])/...
      (sum(Ki([1 2 4 5])./Wi([1 2 4 5]))/(1-Ki(3)));
  end
  av12 = sum(d5spins .* repmat(Ni, NPOINTS, 1), 2);
  
  [wakeProxy,ind1,ind2] = getProxyWake();
  if isempty(wakeProxy), continue, end
  
  flagSWSpin1 = isSW(idxSpin1); flagSWSpin2 = isSW(idxSpin2);
  [wake1,wake2, wakeDescIdx1, d12, d12Plot] = getFinalWake();
  if isempty(wake1) && isempty(wake2), continue, end
  wakeDescIdx1(1) = wakeDescIdx1(1)-expPhase; % probe phase of max wake
  wakedesc(find(idx==idx0,1,'first'), 2:4) = wakeDescIdx1;
  
  wake = wakeModel(idxSpin);
  wExprap = imag(wake); wake = real(wake);
  if ~isempty(wake1)
    if any(abs(wake(idx1))>0), wake = 0.5*(wake+wake1); % average with prev spin
    else, wake = wake1;
    end
    wExprap(idx1) = 0;
  end
  if ~isempty(wake2)
    wake = wake + wake2;
    wExprap(idx2) = 0;
  end
  wakeModel(idxSpin) = wake + 1i*wExprap;
  n_corrected = n_corrected + 1;
  
  % Preliminary apply to the next spin, which can be owerwitten on next step
  wakeModel(idxSpin+NPOINTS) =  real(wake)*1i;
  % Preliminary apply to the spin after the next one, if Solar Wind
  %idx2plus = idxCurrSpin(1) -1 +NPOINTS*2 +idx1;
  %if isSW(idx2plus)
  %  wakeModel(idx2plus) =  real(wake(idx1))*1i;
  %end
  WAMP_THRESHOLD = 1e-4;
  % Apply to prev spin if no wake was found there
  idx1minus = idxSpin1 -NPOINTS;
  if ~any(abs(wakeModel(idx1minus))>WAMP_THRESHOLD) %% chk if 0 is a good threshold;
    wakeModel(idx1minus) = real(wake(idx1))*1i;
  end
  % Apply two spins back, if Solar Wind
  %idx2minus = idxCurrSpin(1) -1 -NPOINTS*2 +idx1;
  %if isSW(idx2minus) && ~any(abs(wakeModel(idx2minus))>WAMP_THRESHOLD)
  %  wakeModel(idx2minus) =  real(wake(idx1))*1i;
  %end
  
  % If no wake1, but ok wake2 from the current and prev spin, copy
  % extrapolated (imaginary) wake1 to next spin
  if ~isempty(wake2) && isempty(wake1) && ...
      any(abs(real(wakeModel(idxSpin2 -NPOINTS)))>WAMP_THRESHOLD) % prev wake2 was good
    wakeModel(idxSpin1 +NPOINTS) = wakeModel(idxSpin1);
  end
  
  plotNow()
  if plotflag_now && idx~=idx0(end)
    plot_step = irf_ask('Step? (0-continue, -1 return) [%]>','plot_step',1);
    if plot_step==0, plotflag = 0;
    elseif plot_step<0, return
    end
  end
  
  % Fist spins of the segment
  if idx2deg(idx)<360
    wakeModel(idxSpin-NPOINTS) = wake;
    wakeModel(idxSpin-NPOINTS*2) = wake;
    n_corrected = n_corrected + 4;
  end
  % Last spins of the segment
  if length(wakeModel)-idx<NPOINTS*(NWSPINS+.5)-1
    wakeModel(idxSpin+NPOINTS) = wake;
    wakeModel(idxSpin+NPOINTS*2) = wake;
    n_corrected = n_corrected + 4;
    break
  end
end


epoch0 = timeIn(1); epochFixedPhaTmp = double(epochFixedPha-epoch0);
epochTmp = double(timeIn-epoch0);
if ~isreal(wakeModel)
  wakeModel = real(wakeModel) + imag(wakeModel);
end
wakeModelOut = interp1(epochFixedPhaTmp,wakeModel,epochTmp,'spline');

irf.log('notice', ['Corrected ', num2str(n_corrected), ' out of ', ...
  num2str(length(idx0)), ' spins.']);
return
% End MAIN

  function res = isSW(idx)
    res = [];
    if isempty(swFlag), return; end
    if any(swFlag(idx)), res = true; else, res = false; end
  end

  function plotNow()
    % Do we need to plot the current interval?
    if plotflag, plot_i = plot_i + 1; end
    if plotflag && plot_i == plot_step
      plotflag_now = true; plot_i = 0;
    else, plotflag_now = false;
    end
    if ~plotflag_now, return, end
    
    h = irf_plot(4,'reset');
    phaPlot = unwrap(pha5spins(:,3)*pi/180)*180/pi - expPhase;
    off = median(d5spins(:,3));
    plot(h(1),phaPlot,d5spins(:,1), 'g',...
      phaPlot,d5spins(:,2), 'g',...
      phaPlot,d5spins(:,4), 'g',...
      phaPlot,d5spins(:,5), 'g',...
      phaPlot,d5spins(:,3), 'b',...
      phaPlot, av12, 'k',...
      phaPlot(ind1)*[1 1], [-2 2]+off, 'm',...
      phaPlot(ind2)*[1 1], [-2 2]+off, 'm',...
      phaPlot, d5spins(:,3)-wake,'r');
    ylabel(h(1),[pair ' [mV/m]']);
    title(h(1),tStUTC)
    set(h(1),'XAxisLocation','top')
    
    plot(h(2),phaPlot,d12Plot,'g',phaPlot, d12,'b');
    ylabel(h(2),['D2(' num2str(pair) ') [mV/m]']);
    
    plot(h(3),phaPlot, wake)
    ylabel(h(3),'Wake [mV/m]');
    set(h(1:3),'XLim',phaPlot([1 end]),'XTick',-360:45:360*2)
    
    wakeModelTmp = zeros(length(idxSpin),5);
    for i=1:5
      ii = idxSpin + (i-3)*NPOINTS;
      wakeModelTmp(:,i) = real(wakeModel(ii)) + imag(wakeModel(ii));
    end
    plot(h(4),phaPlot, d5spins(:,3), 'k', phaPlot, d5spins(:,3)-wake,'r',...
      phaPlot, d5spins(:,3)-wakeModelTmp(:,3),'m');
    ylabel(h(4),[pair ' [mV/m]']);
    hold(h(4), 'on')
    for i=1:5
      if i==3, continue, end
      plot(h(4),phaPlot+360*(i-3), d5spins(:,i), 'k', ...
        phaPlot+360*(i-3), d5spins(:,i)-wakeModelTmp(:,i),'b')
    end
    
    set(h(4),'XLim',phaPlot([1 end])+[-720; 720]), set(h(2:3),'XTickLabel',[])
    set(h(4),'XTick',-360*3:90:360*4)
    set(h(4),'XTickLabel',mod(get(h(4),'XTick'),360))
  end % plotNow

  function [wake,ind1,ind2]  = getProxyWake()
    % First find a proxy wake fit
    % Identify wakes by max second derivative
    
    if isempty(flagSWSpin), FACTOR_SW_AMP = 1; FACTOR_SW_ANG = 1;
    elseif flagSWSpin, FACTOR_SW_AMP = 1.3; FACTOR_SW_ANG = 1.3;
    else, FACTOR_SW_AMP = 0.8; FACTOR_SW_ANG = 0.8;
    end
    
    wake = [];
    d12 = [av12(1)-av12(end); diff(av12)];
    d12 = [d12(1)-d12(end); diff(d12)];
    % Average with 7 points to minimize danger of detecting a wrong maximum
    d12 = w_ave(d12, 7, NPOINTS);
    
    % Ensure wake is symmetrical (with 180 +/-5 deg difference) between the
    % two probes.
    
    if fixedPha(idx) == 0
      ind1 = find(d12 == max(d12(expPhaseIdx))) -1;
      ind2 = find(d12 == min(d12(expPhaseIdx+deg2idx(180)))) -1;
    elseif fixedPha(idx) == 180
      ind1 = find(d12 == min(d12(expPhaseIdx))) -1;
      ind2 = find(d12 == max(d12(expPhaseIdx+deg2idx(180)))) -1;
    end
    if isempty(ind1) || isempty(ind2)
      irf.log('debug', ['No wake at ' tStUTC]);
      return
    end
    if abs(ind2-ind1-deg2idx(180)) > 15*FACTOR_SW_ANG
      irf.log('debug', ['Wake displaced by '...
        num2str(abs(idx2deg(ind2-ind1)-180)') ' deg at ' tStUTC]);
      return
    end
    
    % The proxy wake is naroow (1/2 of the final fit)
    wake_width = deg2idx(WAKE_MAX_HALFWIDTH/2);
    i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
    i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;
    
    % The proxy wake is symmetric
    dav = (d12(i1)-d12(i2))/2;
    cdav = cumsum(dav);
    cdav = cdav - mean(cdav);
    ccdav = cumsum(cdav);
    
    % Wake half-width
    iiTmp = find(abs(ccdav)<max(abs(ccdav))/2);
    wampl = min(iiTmp(iiTmp>23))-max(iiTmp(iiTmp<23)); % XXX: what are these 23s??
    if isempty(wampl) || wampl<WAKE_MIN_AMPLITUDE/FACTOR_SW_AMP
      irf.log('debug', ['Proxy wake is too small at ' tStUTC]);
      return
    end
    
    wake = zeros(NPOINTS,1);
    wake( i1 ) = ccdav;
    wake( i2 ) = -ccdav;
  end % getProxyWake

  function [wake1,wake2,wakedesc1,d12,d12Plot] = getFinalWake()
    % Find final wake shape
    
    % Correct for the proxy wake
    av12_corr = av12 - wakeProxy;
    
    % Find the ground tone and remove it from the data
    x = fft(av12_corr);
    x(3:end-1) = 0;
    av12_corr = av12 -ifft(x, 'symmetric');
    
    % Now find the final fit
    d12 = [av12_corr(1)-av12_corr(end); diff(av12_corr)];
    d12 = [d12(1)-d12(end); diff(d12)];
    d12Plot = d12; % save for plotting
    % Average with only 5 points to get a more fine fit
    d12 = w_ave(d12, 5, NPOINTS);
    
    wake_width = deg2idx(WAKE_MAX_HALFWIDTH);
    i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
    i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;
    
    if fixedPha(idx) == 0, w1 = 'min'; w2 = 'max';
    elseif fixedPha(idx) == 180, w1 = 'max';  w2 = 'min';
    end
    
    [wake1, wakedesc1] = getOneWake(i1,w1,flagSWSpin1);
    wake2 = getOneWake(i2,w2,flagSWSpin2);
    
    function [wake, wakeDesc] = getOneWake(idx,wMode,flag_solar_wind)
      wake = [];
      wakeDesc = NaN(1,3); % to store [phase, amplitude, half width]
      if nargin<2, flag_solar_wind = []; end
      
      if isempty(flag_solar_wind), FACTOR_SW_AMP = 1; FACTOR_SW_SHAPE = 1;
      elseif flag_solar_wind, FACTOR_SW_AMP = 1.3; FACTOR_SW_SHAPE = 1.3;
      else, FACTOR_SW_AMP = .8; FACTOR_SW_SHAPE = 0.8;
      end
      
      % Allow the final fit to be asymmetric
      cdav = cumsum(d12(idx));
      cdav = cdav - mean(cdav);
      ccdav = cumsum(cdav);
      switch wMode
        case 'max', ccdav(ccdav<0) = 0;
        case 'min', ccdav(ccdav>0) = 0;
        otherwise
          error('bad mode for getOneWake')
      end
      ccdav = crop_wake(ccdav);
      wakeDesc(2) = max(abs(ccdav)); % save amplitude of wake
      if wakeDesc(2)< WAKE_MIN_AMPLITUDE/FACTOR_SW_AMP ||...
          wakeDesc(2)>WAKE_MAX_AMPLITUDE
        irf.log('debug', ...
          sprintf('Wake too small/big(%.2f mV/m) at %s', ...
          wakeDesc(2), tStUTC));
        return
      end
      
      if ~isGoodShape(ccdav,FACTOR_SW_SHAPE)
        irf.log('debug', ['Wrong wake shape at ' tStUTC]);
        return
      end
      
      % Wake half-width
      ii =    find( abs(ccdav) <  max(abs(ccdav))/2 );
      iimax = find( abs(ccdav) == max(abs(ccdav))   );
      hw = idx2deg(min(ii(ii>iimax))-max(ii(ii<iimax)));
      if isempty(hw)
        irf.log('debug', ['wrong wake shape at ', ...
          tStUTC ' (spike corner case)']);
        return
      elseif hw < WAKE_MIN_HALFWIDTH/FACTOR_SW_AMP
        irf.log('debug', sprintf('wake is too narrow (%d deg) at %s', ...
          hw, tStUTC));
        return
      end
      wakeDesc(1) = idx2deg(idx(iimax)); % save angle of max wake
      wakeDesc(3) = hw; % save half width
      wake = zeros(NPOINTS,1); wake( idx ) = ccdav;
    end
  end

  function idx = deg2idx(deg)
    idx = round(deg/360*NPOINTS);
  end

  function deg = idx2deg(idx)
    deg = idx/NPOINTS*360;
  end

  function wake = crop_wake(wake)
    % Crop wake side lobes below a defined fraction of the maximum.
    % Use spline interpoltion to reach smooth transition to the zero level
    % outside. If lobes are not located or too close to the beggining or
    % end of the wake segment it is returned unaltered.
    
    %DEBUG=false;
    AMP_FRAC = 0.1; % fraction of amplitude below which we neew to crop
    GAP_WIDTH = round(4*NPOINTS/360); % number of points ower which the wake is required to reach zero
    lenWake = length(wake);
    
    idxx = (1:lenWake)';
    imax = find(abs(wake)==max(abs(wake)));
    wamp = wake(imax);
    
    if wamp<0, iout = wake>wamp*AMP_FRAC;
    else, iout = wake<wamp*AMP_FRAC;
    end
    
    ist = find(((idxx<imax) & iout),1,'last');
    ien = find(((idxx>imax) & iout),1,'first');
    if isempty(ist) || isempty(ien)
      irf.log('debug', 'Avoiding cropping wake, empty "ist" or "ein".');
      return
    elseif (ist<=GAP_WIDTH) || (ien>=lenWake-GAP_WIDTH)
      irf.log('debug', 'Avoiding cropping wake, too close to beginning/end.');
      return
    end
    %if DEBUG
    %  figure;
    %  subplot(2,1,1);
    %  plot(idx, wake, '-black', imax, wamp, '-rO');
    %end
    wake(idxx>=ien+GAP_WIDTH+1 | idxx<=ist-GAP_WIDTH-1) = 0;
    
    iexcl = [ist-GAP_WIDTH:ist, ien:ien+GAP_WIDTH]; % indeces over which to interpolate
    itmp = setxor(idxx,iexcl);
    
    % Pad with zeros at the edges before interpolating
    wakeTmp = interp1([-1; 0; itmp; lenWake+1; lenWake+2],[0; 0; wake(itmp); 0; 0],...
      [-1; 0; idxx; lenWake+1; lenWake+2],'spline');
    wake = wakeTmp(3:end-2);
    
    %if DEBUG
    %  subplot(2,1,1);
    %  hold on;
    %  plot(idx, wake, '-green');
    %  legend('wake original', 'wake max', 'wake edges set to zero');
    %  ylabel('Wake [mV/m]');
    %  subplot(2,1,2);
    %  plot(idx, iout, '-blue*', ist, 1,'-rO', ien, 1, '-rO', itmp, ones(size(itmp)), '-greenO');
    %  legend('iout', 'ist', 'iend', 'itmp');
    %  ylim([-0.1 1.1]); ylabel('indicator true/false');
    %end
  end
end

function av = w_ave(x, np, NPOINTS)
% Weighted average
narginchk(3, 3);
av = zeros(size(x));
if np==7
  m = [.07; 0.15; 0.18; 0.2; 0.18; 0.15; 0.07];
  idx = -3:1:3;
else
  m = [0.1; 0.25; 0.3; 0.25; 0.1];
  idx = -2:1:2;
end
MIDX = max(idx);
for j=1:length(x)
  ii = j + (idx);
  if j<=MIDX
    ii(ii<1) = ii(ii<1) + NPOINTS;
  end
  if j>length(x)-MIDX
    ii(ii>NPOINTS) = ii(ii>NPOINTS) - NPOINTS;
  end
  av(j) = sum(x(ii).*m);
end
end

function res = isGoodShape(s,FACTOR_SW)
% check for shape of the wake fit
RATIO = 0.4; % (Cluster was 0.3)
% ThoNi: one testrun with 20170508 mms1 fast got a spike in frequency in
% the interval 0.3->0.4 compared with intervals 0.4->0.5, 0.5->0.6 etc.
% Therefor try increasing the permitted Ratio to 0.4 compared with 0.3
% which was used for Cluster.
res = true;
if max(s)~=max(abs(s))
  s = -s;
end
maxmax = max(s); % global maxima
d1 = diff(s);
imax = find( (d1(1:end-1).*d1(2:end))<0 ) + 1;
if length(imax)==1
  if s(imax)==maxmax
    return
  else
    % Signal has a minumum
    irf.log('debug', 'Signal has a minumum instead of a maximum.');
    res = false;
    return
  end
elseif isempty(imax)
  % Signal is monotonic
  irf.log('debug', 'Signal is monotonic.');
  res = false;
  return
end
maxima = abs(s(imax));
smax = max( maxima(maxima < maxmax) ); % Second maxima
if smax>maxmax*RATIO*FACTOR_SW
  res = false;
  irf.log('debug', ...
    sprintf('BAD FIT: second max is %0.2f of the main max', smax/maxmax) );
end


end

