function [HBIASSA,wakedesc] = c_efw_hbias_satur(e,pair,phase_2,whip,plotflag)
%C_EFW_HBIAS_SATUR  Check for EFW saturation due to high bias current
%
% [HBIASSA,spikedesc] = c_efw_hbias_satur(e,pair,phase_2,whip,plotflag)
%
% Input:
%   e        - raw EFW data (wE?p12/34)
%   pair     - probe pair (12/34/32)
%   phase_2  - Phase 2
%   plotflag - 0 no plotting of debug
%
% Output:
%   HBIASSA - dirty intervals
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


DEBUG = 0;

HBIASSA = [];
wakedesc=[];

narginchk(3,5)
if nargin<5, plotflag = 0; end
if nargin<4, whip = []; end

if pair~=12 && pair~=32 && pair~=34, error('PAIR must be one of: 12, 32, 34'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .9;
MAX_SPIN_PERIOD = 4.3;
WAKE_MAX_WIDTH = 150; % degrees
WAKE_MIN_WIDTH = 20;  % degrees
WAKE_MAX_CENTER_OFF = 8; % degrees
WAKE_MIN_AMPLITUDE = 5; % mV/m
WAKE_MAX_AMPLITUDE = 780; % mV/m
WAKE_MIN_MAX_RATIO = 5; % Min ration of the max1/max ratio
WAKE_INT_AMPLITUDE = 100; % mV/m
WAKE_MAX_ADC_OFF = 50;
WIDE_WAKE_MIN_AMPLITUDE = 50; % Min amplitude for wide wakes
WIDE_WAKE_MIN_WIDTH = 100;  % Width at which to consider a wake as wide


plot_step = 1;
plot_i = 0;

data_corr = e;
if ~isempty(whip)
  irf_log('proc','blanking Whisper pulses')
  data_corr = caa_rm_blankt(e,whip);
end

i0 = find(phase_2(:,2)==0);
if isempty(i0), irf_log('proc','empty phase'), return, end
if length(i0)<5, irf_log('proc','not enough spins for correction'), return, end

% Hack to add phase=360 at the end and phase=0 at the start of the interval
% in order to try to correct these incomplete spins
if phase_2(end,2)~=360 && ...
    phase_2(i0(end),1)-phase_2(i0(end-1),1)<MAX_SPIN_PERIOD
  phase_2 = [phase_2; phase_2(i0(end),1)*2-phase_2(i0(end-1),1) 0];
end
if phase_2(1,2)~=0 && phase_2(i0(2),1)-phase_2(i0(1),1)<MAX_SPIN_PERIOD
  phase_2 = [phase_2(i0(1),1)*2-phase_2(i0(2),1) 0; phase_2];
  i0 = find(phase_2(:,2)==0);
end

n_spins = length(i0);
NPOINTS = 361;
tt = zeros(NPOINTS,n_spins);
wakedesc = zeros(n_spins, 3)*NaN;
ttime = tt;
sf = c_efw_fsample(e,'hx');
iok = [];

for in = 1:n_spins
  ts = phase_2(i0(in),1);
  i360 = find( phase_2(:,1)>ts & phase_2(:,1)<ts+MAX_SPIN_PERIOD & ...
    phase_2(:,2)==360);
  if isempty(i360)
    irf_log('proc',['gap in phase at ' epoch2iso(ts,1)])
    te = ts + 4.0;
    empty = 1;
  else
    if length(i360)~=1, irf_log('proc',['bogus phase at ' epoch2iso(ts,1)]), end
    te = phase_2(i360(end),1);
    empty = 0;
  end
  
  ttime(:,in) = (ts + (0:1:360) *(te-ts)/360.0)';
  if empty, tt(:,in) = NaN;
  else
    MARG = 0.05;
    eind = find((data_corr(:,1) > ts-MARG) & (data_corr(:,1) < te+MARG));
    eind(isnan(data_corr(eind,2))) = [];
    % Check for data gaps inside one spin.
    if sf>0 && length(eind)<N_EMPTY*(te-ts +MARG*2)*sf
      if DEBUG, irf_log('proc',['data gap at ' epoch2iso(ts,1)]), end %#ok<UNRCH>
      tt(:,in) = NaN;
    else
      if sf==450, dtmp = irf_resamp(data_corr(eind,:), ttime(:,in));
      else
        ddt = 1/sf;
        dtmp = data_corr(eind,:);
        % We linearly extrapolate missing data at at edges
        if dtmp(1,1)>=ts+ddt
          % We miss points at the beginning of the spin
          nm = ceil( (dtmp(1,1)-ts)/ddt );
          t_temp = dtmp(1,1) + ((1:nm) - nm-1)*ddt;
          data_temp = irf_resamp(dtmp,[t_temp'; dtmp(1:2,1)],...
            'linear');
          dtmp = [data_temp(1:end-2,:); dtmp]; %#ok<AGROW>
          clear nm t_temp data_temp
        end
        if dtmp(end,1)<=te-ddt
          % We miss points at the end of the spin
          nm = ceil( (te-dtmp(end,1))/ddt );
          t_temp = dtmp(end,1) + (1:nm)*ddt;
          data_temp = irf_resamp(dtmp,...
            [dtmp(end-1:end,1); t_temp'], 'linear');
          dtmp = [dtmp; data_temp(3:end,:)]; %#ok<AGROW>
          clear nm t_temp data_temp
        end
        dtmp = irf_resamp(dtmp, ttime(:,in), 'spline');
      end
      % Fill small gaps (at edges only?) with zeroes
      % This has a minor influence on the correction procedure
      dtmp(isnan(dtmp(:,2)),2) = 0;
      tt(:,in) = dtmp(:,2);
      clear dtmp
    end
    
    iok = [iok in];
  end
end

clear data_corr

if isempty(iok)
  irf_log('proc','not enough spins for correction')
  return
end

if pair==12, da = 45;
elseif pair==34, da = -45;
end

for in = iok
  % Do we need to plot the current interval?
  if plotflag, plot_i = plot_i + 1; end
  if plotflag && plot_i == plot_step
    plotflag_now = 1;
    plot_i = 0;
  else, plotflag_now = 0;
  end
  
  wakedesc(in,1) = ttime(1,in);
  
  if any(isnan(tt(:,in)))
    if DEBUG, disp(['gap at spin # ' num2str(in)]), end %#ok<UNRCH>
    continue
  end
  
  [min1,imin1,wmin1,max1,imax1,wmax1,min2,max2] = find4minmax(tt(:, in));
  if abs(min1)>abs(max1), w = wmin1; im = imin1;
  else, w = wmax1; im = imax1;
  end
  mm = max(abs(min1),abs(max1));
  
  if plotflag_now
    clf
    ts = ttime(1,in);
    te = ttime(end,in);
    e_tmp = e(e(:,1)>=ts & e(:,1)<te,:);
    plot(ttime(:,in)-ts, tt(:, in), 'b',...
      ttime(imin1,in)*[1 1]-ts, mm*[-1 1], 'r',...
      ttime(imax1,in)*[1 1]-ts, mm*[-1 1], 'r',...
      e_tmp(:,1)-ts,e_tmp(:,2),'g');
    ylabel('E12 [mV/m]');
    
    ticks = (0:8)*45;
    
    % Fix wake position
    set(gca,'XTick',ttime(ticks+1,in)-ts,'XTickLabel',ticks+da-90)
    title(epoch2iso(ts,1))
    xlabel('angle to the Sun [deg]')
  end
  
  if DEBUG
    if isempty(min2) %#ok<UNRCH>
      disp('No secondary peaks')
    else
      disp(sprintf('Max: %.1f mV/m  Ratio: %.1f  Width: %.1f',...
        mm,mm/max(abs(min2),abs(max2)), w))
    end
  end
  
  if (mm > WAKE_MAX_AMPLITUDE) || (mm > WAKE_INT_AMPLITUDE && w >=WAKE_MIN_WIDTH)
    if DEBUG, disp('Matches'), end %#ok<UNRCH>
    wakedesc(in,2) = mm;
    wakedesc(in,3) = w;
  elseif mean(tt(:, in)) > WAKE_MAX_ADC_OFF
    if DEBUG, disp('Possibly matches based on ADC offset'), end %#ok<UNRCH>
    wakedesc(in,2) = -mm;
    wakedesc(in,3) = -w;
  elseif pair~=32	 % Wake position works only for p12 and p34
    if(mm > WAKE_MIN_AMPLITUDE && ~isempty(min2) && ...
        w >= WAKE_MIN_WIDTH)
      % Check for wake position
      im = im +da -90; % Angle with respect to the Sun
      if round( im/180 ) >= 1, im = im -180; end
      if abs(im) > WAKE_MAX_CENTER_OFF
        if DEBUG
          disp(sprintf('wake center offset by : %d degress',im)) %#ok<UNRCH>
        end
      else
        if(w <= WAKE_MAX_WIDTH)
          if( mm/max(abs(min2),abs(max2)) > WAKE_MIN_MAX_RATIO)
            if DEBUG, disp('Matches'), end %#ok<UNRCH>
            wakedesc(in,2) = mm;
            wakedesc(in,3) = w;
          else
            if (w < WIDE_WAKE_MIN_WIDTH) || (mm > WIDE_WAKE_MIN_AMPLITUDE)
              if DEBUG, disp('Possibly matches, but only one spike detected.'), end %#ok<UNRCH>
              wakedesc(in,2) = -mm;
              wakedesc(in,3) = -w;
            end
          end
        end
      end
    else
      if DEBUG, disp(['no saturation at spin # ' num2str(in)]), end %#ok<UNRCH>
    end
  else
    if DEBUG, disp('wake position works only for p12 and p34'), end %#ok<UNRCH>
  end
  if plotflag_now,input('press enter'), end
end

% Take care of any "maybe" matches (i.e. detected one spike only)
maybes=find(wakedesc(iok,2)<0);
if ~isempty(maybes)
  if length(iok) > 1
    neighbor1=maybes-1;
    neighbor2=maybes+1;
    if neighbor1(1) < 1, neighbor1(1)=2; end
    if neighbor2(end) > length(iok), neighbor2(end) = length(iok)-1; end
    idx1=find( isfinite(wakedesc(iok(neighbor1),2)) |  isfinite(wakedesc(iok(neighbor2),2)));
    idx2=find(~isfinite(wakedesc(iok(neighbor1),2)) & ~isfinite(wakedesc(iok(neighbor2),2)));
    if ~isempty(idx1)
      wakedesc(iok(maybes(idx1)),2)= -wakedesc(iok(maybes(idx1)),2);
      wakedesc(iok(maybes(idx1)),3)= -wakedesc(iok(maybes(idx1)),3);
    end
    if ~isempty(idx2), wakedesc(iok(maybes(idx2)),2:3)= NaN; end
  else
    wakedesc(iok,2:3)= NaN;
  end
end

% Join intervals
st = []; et = [];
for in = iok
  if isempty(st)
    if isnan(wakedesc(in,1)) || isnan(wakedesc(in,2)), continue, end
    
    st = wakedesc(in,1) - 4;
    
    % Check if the prev dirty interval is only one spin away
    if ~isempty(HBIASSA) && HBIASSA(end,2) > st -6 % 6 sec is for security
      st = HBIASSA(end,1);
      HBIASSA(end,:) = [];
    else
      % Check for singular points. These must be avoided
      % If next two spins are clean, throw away this one
      if in == iok(end), continue, end
      % Next spin
      if isnan(wakedesc(in+1,2)) && ...
          ( (in == iok(end-1)) || isnan(wakedesc(in+2,2)) )
        continue
      end
    end
  end
  
  if in == iok(end)
    et = wakedesc(in,1) + 8;
    HBIASSA = [HBIASSA; st et];
    continue
  end
  
  if isnan(wakedesc(in+1,2))
    if isempty(et)
      et = wakedesc(in,1) + 4;
    else
      % We are not at the end and the next spin is also clean
      if in+2 == iok(end) || ( in+2 < iok(end) && isnan(wakedesc(in+2,2)))
        HBIASSA = [HBIASSA; st et];
        st = []; et = [];
      else
        et = et + 4;
      end
    end
  else
    et = wakedesc(in,1) + 8;
  end
end

%% help functions
function [min1,imin1,wmin1,max1,imax1,wmax1,min2,max2] = find4minmax(data)

PLOT_FLAG = 0;

if PLOT_FLAG
  clf, plot(data,'k'), hold on, grid on %#ok<UNRCH>
end

data = w_ave(data,5);
if PLOT_FLAG, plot(data); end %#ok<UNRCH>

%d1 = diff(data);
%plot(idx(1:end-1)+.5,d1*2,'r')

min1 = min(data);
max1 = max(data);

imin1 = find(data == min1); imin1 = imin1(1);
imax1 = find(data == max1); imax1 = imax1(1);

if PLOT_FLAG
  plot(imin1,data(imin1),'*',imax1,data(imax1),'*') %#ok<UNRCH>
  %plot(imin1,d1(imin1),'r*',imax1,d1(imax1),'r*')
end

[ii_left, ii_right] = find_peak(data,imin1,0);
wmin1 = ii_right - ii_left +1;

%plot(ii_left,d1(ii_left),'r*',ii_right,d1(ii_right),'r*')

idx_out_min = ii_left:ii_right;
data(idx_out_min) = NaN;
%plot(data,'g')

[ii_left, ii_right] = find_peak(data,imax1,1);
wmax1 = ii_right - ii_left +1;

%plot(ii_left,d1(ii_left),'r*',ii_right,d1(ii_right),'r*')

idx_out_max = ii_left:ii_right;
data(idx_out_max) = NaN;
%plot(data,'m')

d_tmp = data(~isnan(data));

if isempty(d_tmp)
  max2 = [];
  min2 = [];
  %disp('no secondary peaks')
  return
end

max2 = max(d_tmp);
min2 = min(d_tmp);

imin2 = find(data == min2);
imax2 = find(data == max2);

if PLOT_FLAG
  plot(imin2,data(imin2),'o',imax2,data(imax2),'o') %#ok<UNRCH>
  
  disp(sprintf('Max 1: %.1f (%d points)  Max 2: %.1f  Min1/Min2: %.1f',...
    max1,wmax1,max2,max1/max2))
  disp(sprintf('Min 1: %.1f (%d points)  Min 2: %.1f  Max1/Max2: %.1f',...
    min1,wmin1,min2,min1/min2))
end

%% find_peak
function [ii_left, ii_right] = find_peak(data,im,mode)
% mode = 0  minimum
% mode = 1  maximum
N_HALFWIDTHS=2.0;
if ~mode, data = -data; end

ndata = length(data) +1;
N_HALFWIDTHS=N_HALFWIDTHS-1;

ii_left = 1;
ii_right = ndata;
if im > 1
  idxtmp = 1:im-1;
  ii_left = find(data(idxtmp) < 0.5*data(im));
  if isempty(ii_left), ii_left = 1;
  else
    ii_left = idxtmp(max(ii_left))+1;
    ii_left = ii_left-N_HALFWIDTHS*(im-ii_left);
    if ii_left<1, ii_left=1; end
  end
end
if im < ndata
  idxtmp =im:ndata-1;
  ii_right = find(data(idxtmp) < 0.5*data(im));
  if isempty(ii_right), ii_right = ndata;
  else
    ii_right = idxtmp(min(ii_right)) -1;
    ii_right = ii_right+N_HALFWIDTHS*(ii_right-im);
    if ii_right > (ndata-1), ii_right=ndata-1; end
  end
end

%% w_ave
function av = w_ave(x,np)
% Weighted average
NPOINTS = length(x);
if nargin<2, np=5; end
av = zeros(size(x));
if np==7
  m = [.07 .15 .18 .2 .18 .15 .07]';
  idx = -3:1:3;
else
  m = [.1 .25 .3 .25 .1]';
  idx = -2:1:2;
end
MIDX = max(idx);

for j=1:length(x)
  ii = j + (idx);
  if j<=MIDX, ii(ii<1) = ii(ii<1) +NPOINTS; end
  if j>length(x)-MIDX, ii(ii>NPOINTS) = ii(ii>NPOINTS) -NPOINTS; end
  av(j) = sum(x(ii).*m);
end
