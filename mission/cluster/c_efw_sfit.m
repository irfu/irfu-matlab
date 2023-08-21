function spinfit = c_efw_sfit(pair,fout,maxit,minpts,te,data,tp,ph,method,tmmode)
%C_EFW_SFIT  Produce spin fit data values (ex, ey)
%   of EFW data frome given probe pair at 4 second resolution.
%
% SPINFIT = C_EFW_SFIT(PAIR,FOUT,MAXIT,MINPTS,TE,DATA,TP,PH,METHOD,[TMMODE])
%
% Input:
%  PAIR - probe pair used (12, 32, 34)
%  FOUT - minimum fraction of fit standard deviation that defines an outlier
%         (zero means no removal of outliers). Has no effect if METHOD=1
%  MAXIT - maximum number of iterations (zero means infinity)
%  MINPTS - minimum number of data points to perform fit
%      (set to 5 if smaller number if given)
%  TE - EFW time in seconds (isGetDataLite time)
%  DATA - EFW data from pair in mV/m, should correspond to te
%  TP - Ephemeris time in seconds (isGetDataLite time)
%  PH - Ephemeris phase_2 in degr (sun angle for s/c Y axis), should
%      correspond to tp
%  METHOD - 0: c_efw_onesfit, Matlab routine by AIE
%           1: c_efw_spinfit_mx (default), BHN Fortran source obtained from KTH
%  TMMODE - 'hx' (default) or 'ib', used for determination of sampling freq
%
% Output:
%  SPINFIT = [TS,EX,EY,OFFSET,SDEV0,SDEV,ITER,NOUT]
%  TS - time vector in seconds
%  EX - E-field x-component in DSI coordinates (almost GSE)
%  EY - E-field y-component in DSI coordinates (almost GSE)
%  OFFSET - mean value of input data
%  SDEV0 - standard deviation in first fit. Has no meaning if METHOD=1
%  SDEV - standard deviation in final fit
%  ITER - number of iterations (one if OK at once)
%  NOUT - number of outliers removed
%
% This function chops up the time series in four second
% intervals, each of which are analysed by c_efw_onesfit.
% Spins with less then 90% of data are disregarded.
%
% Example: (assuming valid db exist, see isGetDataLite)
% To plot spinfits of data from probe pair 34 (obtained during
%    the first 10 minutes of the Donald Duck show on Swedish
%    national TV on Christmas Eve 2002), defining outliers as
%    points outside 3 standard deviations, allowing maximum 10
%    iterations for each spin, and demanding 20 data points for
%    the fit to be valid, do as follows:
%  [t34,e34] = isGetDataLite(db, [2002 12 24 14 00 00], 600, ...
%  'Cluster', '4', 'efw','E','p34','10Hz','hx');
%  [tpha,pha] = isGetDataLite(db, [2002 12 24 14 00 00], 600, ...
%  'Cluster', '4', 'ephemeris','phase_2','','','');
%  sp34 = c_efw_sfit(34,3,10,20,t34,double(e34),tpha,double(pha));
%  t0 = toepoch([2002 12 24 14 00 00]);
%  t = sp34(:,1) - t0;
%  ex = sp34(:,2);
%  ey = sp34(:,3);
%  plot(t,ex,'k',t,ey,'r');
%  xlabel('Time from 2002-12-24 14:00:00 [s]');
%  ylabel('E [mV/m]');
%  title('Cluster SC4 EFW Spin Fits (Ex black, Ey red)')
%
% See also C_PHASE, C_EFW_ONESFIT

narginchk(8,10)

if pair~=12 && pair~=32 && pair~=34 && pair~=42, error('PAIR must be one of: 12, 34, 32, 42'), end

if ~isequal(size(te),size(data))
  error('TE and DATA vectors must be the same size.')
end
if ~isequal(size(tp),size(ph))
  error('TP and PH vectors must be the same size.')
end

if nargin<10
  tmmode='hx';
end

% Check if the phase is good
[ie,ip] = irf_find_comm_idx(te,tp);
if isempty(ie) || isempty(ip)
  % Try to compute phase
  irf_log('proc','computing phase with c_phase')
  aa = c_phase(te,[tp ph]);
  [ie,ip] = irf_find_comm_idx(te,aa(:,1));
  if isempty(ie) || isempty(ip)
    irf_log('proc','phase/E is empty')
    spinfit = [];
    return
  else
    tp = aa(:,1);
    ph = aa(:,2);
  end
end
te = te(ie); data = data(ie);
tp = tp(ip); ph = ph(ip);
clear ie ip

% Set default method to BHN
if nargin < 9, method = 1; end

if method==1 || method==2
  if exist('c_efw_spinfit_mx','file')~=3
    method = 0;
    disp('cannot find mex file, defaulting to Matlab code.')
  end
end

% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .9;
MAX_SPIN_PERIOD = 4.3;
MIN_SPIN_PERIOD = 3.7;

% Guess the sampling frequency
sf = c_efw_fsample(te,tmmode);

if sf == 0, error('cannot guess the sampling frequency'), end

tpha = tocolumn(tp);
pha = tocolumn(ph);
% Calcluate phase (in rad) at EFW sample times:
pha = unwrap(pi*pha/180);
% Find phase of given pair:
if pair == 12, pha = pha + 3*pi/4;
elseif pair == 32, pha = pha + pi/2;
elseif pair == 34, pha = pha + pi/4;
elseif pair == 42, pha = pha + pi;
else, error('probe pair must be one of 12, 34, 32 or 42')
end


% Do it:
if method==1 || method==2
  te = torow(te);
  data = torow(data);
  ind = find(~isnan(data));
  % Check if we have at least one spin of data
  if length(ind)<N_EMPTY*4*sf
    irf_log('fcal','not enough data points for spinfit');
    spinfit = []; return
  end
  if method==2
    NCOLS = 10;
    NTERM = 5;
  else
    NCOLS = 8;
    NTERM = 3;
  end
  pha = torow(pha(ind));
  [ts,sfit,sdev,iter,nout] = ...
    c_efw_spinfit_mx(maxit,N_EMPTY*4*sf,NTERM,...
    te(ind),data(ind),pha);
  ind = find( sdev~=-159e7 );
  n_gap = length(sdev) -length(ind);
  n = length(sdev);
  spinfit = zeros(n,NCOLS);
  spinfit(:,1) = ts;		        % time
  spinfit(:,2:end) = NaN;
  if ~isempty(ind)
    spinfit(ind,2) = sfit(2,ind);	% Ex
    spinfit(ind,3) = -sfit(3,ind);% Ey, - Because s/c is spinning upside down
    spinfit(ind,4) = sfit(1,ind);
    spinfit(ind,5) = sdev(ind);
    spinfit(ind,6) = sdev(ind);
    spinfit(ind,7) = iter(ind);
    spinfit(ind,8) = nout(ind);
    if method==2
      spinfit(ind,9) = sfit(4,ind); % 2 omega signals
      spinfit(ind,10) = sfit(5,ind);
    end
  end
else
  % Original version using polyfit() by Anders.Eriksson@irfu.se, 13 December 2002

  % Turn off warnings for badly conditioned polynomial:
  warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

  % Chop up time interval
  % We always start at 0,4,8.. secs, so that we have
  % the same timelines an all SC at 2,6,10... sec
  tstart = fix(min(te)/4)*4;
  tend = max(te);

  n = floor((tend - tstart)/4) + 1;
  spinfit = zeros(n,8);

  n_gap = 0;
  for i=1:n
    tsfit = tstart + (i-1)*4 +2;
    ind = find( ( te >= tsfit-2 ) & ( te < tsfit+2 ) );

    % Check for data gaps inside one spin.
    if sf>0 && length(ind)<N_EMPTY*MAX_SPIN_PERIOD*sf
      n_gap = n_gap + 1;
      continue
    end

    % Compute spin period
    pol = polyfit(tpha(ind),pha(ind),1);
    spinp = 2*pi/pol(1);
    if spinp > MAX_SPIN_PERIOD || spinp < MIN_SPIN_PERIOD
      irf_log('proc',sprintf('bad spin period %.4f s',spinp));
      n_gap = n_gap + 1;
      continue
    end

    % Clear NaNs
    ind(isnan(data(ind))) = [];

    % Check for data gaps inside one spin.
    if sf>0 && length(ind)<N_EMPTY*sf*4, continue, end

    % Use Matlab version by AIE
    spinfit(i - n_gap,:) = c_efw_onesfit(pair,fout,maxit,minpts,te(ind),...
      data(ind),te(ind),ph(ind));
    spinfit(i - n_gap, 1) = tsfit;
  end
  spinfit = spinfit(1:n - n_gap, :);

  warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
end

irf_log('proc',sprintf('%d spins processed, %d gaps found',n,n_gap))
