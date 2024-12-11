function [data, n_corrected, wakedesc] = mms_sdp_swwake(e, pair, phase_2, timeIn, sampleRate)
%C_EFW_SWWAKE  Correct raw EFW E data for wake in the solar wind
%
% [data, n_corrected, wakedesc] = mms_sdp_swwake(e, pair, phase_2, timeIn, sampleRate)
%
% Input:
%   e          - raw E-field data
%   pair       - probe pair ('e12' or 'e34')
%   phase_2    - spinphase corresponding to measurements
%   timeIn     - time of measurement (int64, tt2000 ns)
%   sampleRate - nominal sample rate of data (Hz)
%
% Output:
%   data - corrected data
%   n_corrected - numer of spins corrected for wake
%   wakedesc - wake description
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
% one spin takes about 19-20 second (approx phase resolution = 360/(19×32)
% = 0.59 deg per sample).

narginchk(5,5)
global MMS_CONST;
if(isempty(MMS_CONST)), MMS_CONST = mms_constants; end

n_corrected = 0;
data = e;
data_corr = e;
wakedesc = [];

% N_EMPTY 0.75 means that we use only spins with more then 75% points.
N_EMPTY = 0.9; % (Cluster was 0.9)
%MAX_SPIN_PERIOD = 4.3; % sec for Cluster
MAX_SPIN_PERIOD = 10^9*60/MMS_CONST.Spinrate.min; % = 20 sec
WAKE_MAX_HALFWIDTH = 35; % degrees, (Cluster was 45 deg)
WAKE_MIN_HALFWIDTH = 9; %11;  % degrees, (Cluster was 11 deg)
WAKE_MIN_AMPLITUDE = 0.4; % mV/m, (Cluster was 0.4 mV/m)
WAKE_MAX_AMPLITUDE = 7; % mV/m, (Cluster was 7 mV/m)
plot_step = 1;
plot_i = 0;
plotflag = false;

switch pair
  case {'e12', 'e34'}
    % Fix wake position
    expPhase = round((-18:17) + 180/pi*MMS_CONST.Phaseshift.(pair))+1; % Symmetric test
    % As we process each spin separatly (as identified by shift at 360->0)
    % check to see if any expected wake is split between two spins.
    % If so then artificially shift the "phase_2" and "expPhase" by the
    % overlapping amount to ensure it does not occur. This should hopefully
    % avoid having problem with wake removed on a probe pair at spin "N"
    % while not removed from spin "N+1" in what is essentially the same
    % wake on the same probe pair.
    if ( min(expPhase) - WAKE_MAX_HALFWIDTH < 0 )
      DELTA = -(min(expPhase)-WAKE_MAX_HALFWIDTH);
      phase_2 = mod(phase_2 + DELTA, 360);
      expPhase = expPhase + DELTA;
    elseif ( max(expPhase) + 180 + WAKE_MAX_HALFWIDTH > 360 )
      DELTA = -(max(expPhase)+WAKE_MAX_HALFWIDTH-180);
      phase_2 = mod(phase_2 + DELTA, 360);
      expPhase = expPhase + DELTA;
    end
  otherwise
    errStr = 'Pair must be one of: "e12" or "e34"';
    irf.log('critical', errStr);
    error(errStr);
end

% Convert time from ttns (int64) to double keeping for interp1 to work,
% while keeping original input variable "timeIn" (used debug/log messages).
time = double(timeIn-timeIn(1));
epoch0 = EpochTT(timeIn(1)).epochUnix;

% Find numer of spins using the fact that Z-phase is monotonically
% increasing, except for it being modulo 360.
i0 = find(diff(phase_2)<0);
if length(i0)<=5
  irf.log('notice', 'Not enough spins for correction.');
  return
end

% Hack to add phase=360 at the end and phase=0 at the start of the interval
% in order to try to correct these incomplete spins
% if phase_2(end)~=0 && time(i0(end)) - time(i0(end-1)) < MAX_SPIN_PERIOD
%   time = [time; time(i0(end))*2 - time(i0(end-1))];
%   phase_2 = [phase_2; 0];
% end
% if phase_2(1)~=0 && time(i0(2)) - time(i0(1)) < MAX_SPIN_PERIOD
%   time = [time(i0(1))*2 - time(i0(2)); time];
%   phase_2 = [0; phase_2];
%   i0 = find(diff(phase_2)<0);
% end

n_spins = length(i0);
NPOINTS = 361;
tt = zeros(NPOINTS, n_spins);
wakedesc = NaN(n_spins*2, 4);
ttime = tt;
iok = [];

ddt = 10^9/sampleRate; % timeIn is in ns (tt2000), sampleRate in Hz.
MARG = 0.05*10^9;

for in = 1:n_spins
  ts = time(i0(in));
  i360 = find( time > ts & ...
    time < ts+MAX_SPIN_PERIOD & ...
    [diff(phase_2) <= 0; 1]);
  if isempty(i360)
    irf.log('debug', ['gap in phase at ', irf_time(timeIn(i0(in)),'ttns>utc')]);
    te = ts + 20.0; %FIXME + 20 (was "+4" for Cluster)?
    empty = 1;
  else
    if length(i360)~=1
      irf.log('debug', ['bogus phase at ', irf_time(timeIn(i0(in)),'ttns>utc')]);
    end
    te = time(i360(end));
    empty = 0;
  end

  ttime(:, in) = (ts + (0:1:360) *(te-ts)/360.0)';
  if empty
    tt(:, in) = NaN;
  else
    eind = find((time > ts-MARG) & (time < te+MARG));
    eind(isnan(data_corr(eind))) = [];
    % Check for data gaps inside one spin.
    if sampleRate>0 && length(eind)<N_EMPTY*(te-ts +MARG*2)*sampleRate/10^9
      irf.log('debug',['Data gap at ', irf_time(timeIn(i0(in)),'ttns>utc')]);
      tt(:, in) = NaN;
    else
      dtmp = [time(eind), data_corr(eind,:)];
      % We linearly extrapolate missing data at at edges
      if dtmp(1,1)>=ts+ddt
        % We miss points at the beginning of the spin
        nm = ceil( (dtmp(1,1)-ts)/ddt );
        t_temp = dtmp(1,1) + ((1:nm) - nm-1)*ddt;
        data_temp = irf_resamp(dtmp, [t_temp'; dtmp(1:2,1)], 'linear');
        dtmp = [data_temp(1:end-2,:); dtmp]; %#ok<AGROW>
        clear nm t_temp data_temp
      end
      if dtmp(end,1)<=te-ddt
        % We miss points at the end of the spin
        nm = ceil( (te-dtmp(end,1))/ddt );
        t_temp = dtmp(end,1) + (1:nm)*ddt;
        data_temp = irf_resamp(dtmp, [t_temp'; dtmp(end-1:end, 1)], 'linear');
        dtmp = [dtmp; data_temp(3:end,:)]; %#ok<AGROW>
        clear nm t_temp data_temp
      end
      dtmp = irf_resamp(dtmp, ttime(:,in), 'spline');
      % Fill small gaps (at edges only?) with zeroes
      % This has a minor influence on the correction procedure
      dtmp(isnan(dtmp(:,2)), 2) = 0;
      tt(:,in) = dtmp(:,2);
      clear dtmp
    end

    % Identify spins for which we attempt to correct wake
    if (in>=6 && sum(isnan(tt(1,in - (0:1:5) )))<=1)
      % We allow max one data gap within 6 spins
      iok = [iok in-3]; %#ok<AGROW>
    end
    if (in==n_spins || n_spins==5) && ~any(isnan(tt(1,in - (0:1:4) )))
      % Special case when we have only 5 spins or
      % it is the last spin we are working with
      iok = [iok in-2]; %#ok<AGROW>
    end
  end
end

clear data_corr

if isempty(iok)
  irf.log('notice', 'Not enough spins for correction.');
  return
end

prevSpinGood = false; currentSpinGood = false;
for in = iok
  irf.log('debug', sprintf('processing spin nr: %i', in));
  %   if(in == 5), keyboard; irf.log('debug'); end
  prevPrevSpinGood = prevSpinGood;
  prevSpinGood = currentSpinGood; currentSpinGood = true;

  % Do we need to plot the current interval?
  if plotflag
    plot_i = plot_i + 1; %#ok<UNRCH>
  end
  if plotflag && plot_i == plot_step
    plotflag_now = 1; %#ok<UNRCH>
    plot_i = 0;
  else
    plotflag_now = 0;
  end
  ts = ttime(1, in);
  if in==n_spins-2
    % Last interval
    idx = -2:1:2;
  else
    % Check for a data gap
    nans = isnan(tt(1,in + (-2:1:3) ));
    if any(nans)
      idx = 1:6;
      idx = idx(xor(idx,nans)) -3;
    else
      idx = -2:1:2;
    end
  end

  % Spin in the middle has maximum weigth
  av12 = sum(tt(:, in+idx) .* repmat([.1, .25, .3, .25, .1], NPOINTS, 1), 2);

  % First find a proxy wake fit
  % Identify wakes by max second derivative
  d12 = [av12(1)-av12(end); diff(av12)];
  d12 = [d12(1)-d12(end); diff(d12)];
  % Average with 7 points to minimize danger of detecting a wrong maximum
  d12 = w_ave(d12, 7, NPOINTS);

  % Ensure wake is symmetrical (with 180 +/-5 deg difference) between the
  % two probes.
  % ThoNi: one testrun with 20170508 mms1 fast got a majority of displaced
  % wakes at 8 deg (when running with limit +/-5 deg) at times with clearly
  % visible wakes (ie 10:00 UTC). Therefor try to run with +/-8 deg instead
  % of +/-5 deg as used by Cluster.
  % ThoNi: one testrun with 20171115 mms2 fast, allow up to +/-15 deg.
  ind1 = find(d12 == max(d12(expPhase))) -1;
  ind2 = find(d12 == min(d12(expPhase+180))) -1;
  if abs(ind2-ind1-180)>15
    irf.log('debug', ['Wake displaced by ', num2str(abs(ind2-ind1-180)'), ...
      ' deg at ', irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
    wakedesc([in*2-1 in*2], :) = NaN;
    currentSpinGood = false;
  end

  if currentSpinGood
    % The proxy wake is naroow (1/2 of the final fit)
    wake_width = fix(WAKE_MAX_HALFWIDTH/2);
    i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
    i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;

    % The proxy wake is symmetric
    dav = (d12(i1)-d12(i2))/2;
    cdav = cumsum(dav);
    cdav = cdav - mean(cdav);
    ccdav = cumsum(cdav);

    % Save wake description
    fw = (mod(ind2,NPOINTS)+1<mod(ind1,NPOINTS)+1);
    wakedesc(in*2-1+fw, 1)     = ttime(mod(ind1,NPOINTS)+1, in);
    wakedesc(in*2-1+fw, 2)     = ind1;
    wakedesc(in*2-fw, 1)       = ttime(mod(ind2,NPOINTS)+1, in);
    wakedesc(in*2-fw, 2)       = ind2;
    wakedesc([in*2-1 in*2], 3) = max(abs(ccdav));
    % Wake half-width
    ii = find(abs(ccdav)<max(abs(ccdav))/2);
    wampl = min(ii(ii>23))-max(ii(ii<23));
  end

  if currentSpinGood && (isempty(wampl) || wampl<WAKE_MIN_AMPLITUDE)
    irf.log('debug', ['Proxy wake is too small at ', ...
      irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
    wakedesc([in*2-1 in*2], :) = NaN;
    currentSpinGood = false;
  end

  if currentSpinGood
    wakedesc([in*2-1 in*2], 4) = wampl;

    wakeProxy = zeros(NPOINTS,1);
    wakeProxy( i1 ) = ccdav;
    wakeProxy( i2 ) = -ccdav;

    % Correct for the proxy wake
    av12_corr = av12 - wakeProxy;
    % Find the ground tone and remove it from the data
    x = fft(av12_corr);
    x(3:359) = 0;
    av12_corr = av12 -ifft(x, 'symmetric');

    % Now find the final fit
    d12 = [av12_corr(1)-av12_corr(end); diff(av12_corr)];
    d12 = [d12(1)-d12(end); diff(d12)];
    if plotflag_now
      d12_tmp = d12;  %#ok<UNRCH> % save for plotting
    end
    % Average with only 5 points to get a more fine fit
    d12 = w_ave(d12, 5, NPOINTS);

    wake_width = WAKE_MAX_HALFWIDTH;
    i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
    i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;

    % Allow the final fit to be asymmetric
    cdav = cumsum(d12(i1));
    cdav = cdav - mean(cdav);
    ccdav1 = cumsum(cdav);
    ccdav1 = crop_wake(ccdav1);
    cdav = cumsum(d12(i2));
    cdav = cdav - mean(cdav);
    ccdav2 = cumsum(cdav);
    ccdav2 = crop_wake(ccdav2);
  end

  if currentSpinGood && ...
      ( max(max(abs(ccdav1)),max(abs(ccdav2)))< WAKE_MIN_AMPLITUDE ||...
      max(max(abs(ccdav1)),max(abs(ccdav2)))>WAKE_MAX_AMPLITUDE )
    irf.log('debug', ...
      sprintf('Wake too small/big(%.2f mV/m) at %s', ...
      max(max(abs(ccdav1)), max(abs(ccdav2))), ...
      irf_time(int64(ts)+timeIn(1),'ttns>utc')));
    wakedesc([in*2-1 in*2], :) = NaN;
    currentSpinGood = false;
  end
  if currentSpinGood
    [goodShapeCcdav1, maxmax1, smax1] = isGoodShape(ccdav1);
    [goodShapeCcdav2, maxmax2, smax2] = isGoodShape(ccdav2);
    if ~(goodShapeCcdav1 && goodShapeCcdav2)
      irf.log('debug', ['Wrong wake shape at ', ...
        irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
      wakedesc([in*2-1 in*2], :) = NaN;
      currentSpinGood = false;
    end
  end

  if currentSpinGood
    % Save wake description
    wakedesc(in*2-1+fw, 3) = max(abs(ccdav1));
    % Wake half-width
    ii =    find( abs(ccdav1) <  max(abs(ccdav1))/2 );
    iimax = find( abs(ccdav1) == max(abs(ccdav1))   );
    iimax = min(ii(ii>iimax))-max(ii(ii<iimax));
    if ~isempty(iimax)
      wakedesc(in*2-1+fw, 4) = iimax;
    else
      irf.log('debug', ['wrong wake shape at ', ...
        irf_time(int64(ts)+timeIn(1),'ttns>utc'), ' (spike corner case)']);
      wakedesc([in*2-1 in*2], :) = NaN;
      currentSpinGood = false;
      continue
    end
    wakedesc(in*2-fw,3) = max(abs(ccdav2));
    % Wake half-width
    ii =    find( abs(ccdav2) <  max(abs(ccdav2))/2 );
    iimax = find( abs(ccdav2) == max(abs(ccdav2))   );
    iiDiff = min(ii(ii>iimax)) - max(ii(ii<iimax));
    if ~isempty(iiDiff)
      wakedesc(in*2-fw,4) = iiDiff;
    else
      irf.log('debug',['wrong wake shape at ', ...
        irf_time(int64(ts)+timeIn(1),'ttns>utc'), ' (spike corner case)']);
      currentSpinGood = false;
      continue
    end
    clear ii iimax
  end

  if currentSpinGood && ...
      min(wakedesc(in*2-fw,4),wakedesc(in*2-1+fw,4))< WAKE_MIN_HALFWIDTH
    irf.log('debug', sprintf('wake is too narrow (%d deg) at %s', ...
      min(wakedesc(in*2-fw,4), wakedesc(in*2-1+fw,4)), ...
      irf_time(int64(ts)+timeIn(1),'ttns>utc')));
    wakedesc([in*2-1 in*2], :) = NaN;
    currentSpinGood = false;
  end

  if currentSpinGood
    wake = zeros(NPOINTS,1); wake( i1 ) = ccdav1; wake( i2 ) = ccdav2;
    wakePrev = wake;
  elseif prevSpinGood && prevPrevSpinGood
    wake = wakePrev;
    irf.log('debug','using wake shape from the previous spin');
    currentSpinGood = false;
  else, continue
  end

  if plotflag_now, h = irf_plot(4,'reset'); PlotPanels13, end  %#ok<UNRCH>

  % Correct the spin in the middle
  ind = find(time>=ttime(1,in) & time<ttime(end,in));
  if ~isempty(ind)
    wake_e = irf_resamp([ttime(:,in) wake], time(ind));
    data(ind) = data(ind) - wake_e(:,2);
    n_corrected = n_corrected + 1;

    if plotflag_now
      PlotED(h(4)) %#ok<UNRCH>
    end
  end

  % Correct edge spinsat start
  cox = [];
  if in==iok(1) || (in~=iok(1) && in-1~=iok(find(iok==in)-1))
    % If the previous spin was not corrected
    % we correct it here
    if in==iok(1) && in>3
      % We try to correct the spin before the start of the
      % interval which has a data gap
      cox = - (1:3);
      n_corrected = n_corrected + 3;
    else
      cox = - (1:2);
      n_corrected = n_corrected + 2;
    end
    prevSpinGood = true;
  elseif currentSpinGood && prevSpinGood && ~prevPrevSpinGood % Correct previous spin if in was bad
    irf.log('debug','correcting the prev-previous spin')
    cox = -2; n_corrected = n_corrected + 1; wake = wakePrev;
  end

  if in==iok(end) || (in~=iok(end) && in+1~=iok(find(iok==in)+1))
    if ~isempty(cox)
      irf.log('notice', ['single at ', irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
    else
      irf.log('notice', ['stop at ', irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
    end
    if n_spins-in>=3
      % We try to correct the spin at the end of the entire
      % intrval which usually has a data gap
      cox = [cox 1:3]; %#ok<AGROW>
      n_corrected = n_corrected + 3;
    else
      cox = [cox 1:2]; %#ok<AGROW>
      n_corrected = n_corrected + 2;
    end
  elseif ~isempty(cox)
    irf.log('notice', ['start at ', irf_time(int64(ts)+timeIn(1),'ttns>utc')]);
  end
  for cx = cox
    ind = find(time>=ttime(1,in+cx) & time<ttime(end,in+cx));
    if ~isempty(ind)
      irf.log('debug', ['Correcting spin: ', num2str(cx)]);
      wake_e = irf_resamp([ttime(:,in+cx) wake], time(ind));
      data(ind) = data(ind) - wake_e(:,2);
      if plotflag_now
        if ttime(1,in+cx)<ts, ts = ttime(1,in+cx); end
        if ttime(end,in+cx)>te, te = ttime(end,in+cx); end
        hold(h(4),'on')
        PlotED(h(4))
        hold(h(4),'off')
      end
    end
  end
  if plotflag_now && in~=iok(end)
    title(['Plot created: ', char(datetime("now","Format","uuuu/MM/dd")), '. Wakes on pair: ', pair, ' spin: ', num2str(in)]); %#ok<UNRCH>
    plot_step = irf_ask('Step? (0-continue, -1 return) [%]>','plot_step',1);
    if plot_step==0
      plotflag = 0;
    elseif plot_step<0
      return
    end
  end
end

% Save wake position only inside 0-180 degrees
%wakedesc(wakedesc(:,2)>180, 2) = wakedesc(wakedesc(:,2)>180, 2)-180;
wakedesc(isnan(wakedesc(:,1)), :) = [];
% store phase
wakedesc(:,2)=wakedesc(:,2)-expPhase(floor(length(expPhase)/2));

irf.log('notice', ['Corrected ', num2str(n_corrected), ' out of ', ...
  num2str(n_spins), ' spins.']);

  function PlotPanels13
    ts = ttime(1,in);
    te = ttime(end,in);
    plot(h(1),ttime(:,in)-ts, tt(:, in), 'b',...
      ttime(:,in)-ts, tt(:, in + ([-2 -1 1 2]) ), 'g',...
      ttime(:,in)-ts, av12, 'k',...
      ttime(ind1,in)*[1 1]-ts, [-2 2], 'r',...
      ttime(ind2,in)*[1 1]-ts, [-2 2], 'r',...
      ttime(:,in)-ts, av12-wake,'r');
    ylabel(h(1),'E12 [mV/m]');
    irf_timeaxis(h(1),ts); xlabel(h(1),'');
    set(h(1),'XLim',[0 te-ts])

    plot(h(2),ttime(:,in)-ts,d12_tmp,'g',ttime(:,in)-ts, d12,'b');
    ylabel(h(2),['D2(E' num2str(pair) ') [mV/m]']);
    irf_timeaxis(h(2),ts); xlabel(h(2),'');
    set(h(2),'XLim',[0 te-ts])

    plot(h(3),ttime(:,in)-ts, wake, ...
      [ttime(1,in) ttime(end,in)] -ts, [maxmax1 maxmax1],'--g', ...
      [ttime(1,in) ttime(end,in)] -ts, [maxmax2 maxmax2],'--b', ...
      [ttime(1,in) ttime(end,in)] -ts, [smax1 smax1],'--r', ...
      [ttime(1,in) ttime(end,in)] -ts, [smax2 smax2],'--y')
    if strcmp(pair,'e12'), location='west'; else, location='east'; end
    legend(h(3), 'Wake', 'W1 G Max', 'W2 G Max', 'W1 2nd Max', 'W2 2nd Max','location',location);
    ylabel(h(3),'Wake [mV/m]');
    irf_timeaxis(h(3),ts);
    set(h(3),'XLim',[0 te-ts])
  end

  function PlotED(hca)
    irf_plot(hca,{[getEpoch(time(ind)), double(e(ind,:))],[getEpoch(time(ind)), double(data(ind,:))]},'comp')
    ylabel(hca,['E' num2str(pair) ' [mV/m]']);
    irf_zoom(hca,'x',getEpoch([ts te]))
    function epo = getEpoch(t)
      epo = double(t)*1e-9 + epoch0;
    end
  end
end % main

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

function [res, maxmax, smax] = isGoodShape(s)
% check for shape of the wake fit
RATIO = 0.4; % (Cluster was 0.3)
% ThoNi: one testrun with 20170508 mms1 fast got a spike in frequency in
% the interval 0.3->0.4 compared with intervals 0.4->0.5, 0.5->0.6 etc.
% Therefor try increasing the permitted Ratio to 0.4 compared with 0.3
% which was used for Cluster.
res = true; smax=0; negMax=false;
if max(s)~=max(abs(s))
  s = -s;
  negMax=true;
end
maxmax = max(s); % global maxima
d1 = diff(s);
imax = find( (d1(1:end-1).*d1(2:end))<0 ) + 1;
if isscalar(imax)
  if s(imax)==maxmax
    if negMax, maxmax=-maxmax; smax=-smax; end
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
if smax>maxmax*RATIO
  res = false;
  irf.log('debug', ...
    sprintf('BAD FIT: second max is %0.2f of the main max', smax/maxmax) );
end
if negMax
  maxmax=-maxmax; smax=-smax;
end
end

function wake = crop_wake(wake)
% Crop wake side lobes below a defined fraction of the maximum.
% Use spline interpoltion to reach smooth transition to the zero level
% outside. If lobes are not located or too close to the beggining or
% end of the wake segment it is returned unaltered.

%DEBUG=false;
AMP_FRAC = 0.1; % fraction of amplitude bewlo which we neew to crop
GAP_WIDTH = 8; % number of points ower which the wake is required to reach zero
lenWake = length(wake);

idx = (1:lenWake)';
imax = find(abs(wake)==max(abs(wake)));
wamp = wake(imax);

if wamp<0, iout = wake>wamp*AMP_FRAC;
else, iout = wake<wamp*AMP_FRAC;
end

ist = find(((idx<imax) & iout),1,'last');
ien = find(((idx>imax) & iout),1,'first');
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
wake(idx>=ien+GAP_WIDTH+1 | idx<=ist-GAP_WIDTH-1) = 0;

iexcl = [ist-GAP_WIDTH:ist, ien:ien+GAP_WIDTH]; % indeces over which to interpolate
itmp = setxor(idx,iexcl);

% Pad with zeros at the edges before interpolating
wakeTmp = interp1([-1; 0; itmp; lenWake+1; lenWake+2],[0; 0; wake(itmp); 0; 0],...
  [-1; 0; idx; lenWake+1; lenWake+2],'spline');
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
