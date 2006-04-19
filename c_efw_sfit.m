function spinfit = c_efw_sfit(pair,fout,maxit,minpts,te,data,tp,ph,method)
%C_EFW_SFIT produce spin fit data values (ex, ey)
%   of EFW data frome given probe pair at 4 second resolution.
%
% spinfit = c_efw_sfit(pair,fout,maxit,minpts,te,data,tp,ph,method)
%
% Input:
%  pair - probe pair used (12, 32, 34)
%  fout - minimum fraction of fit standard deviation that defines an outlier 
%         (zero means no removal of outliers). Has no effect if METHOD=1
%  maxit - maximum number of iterations (zero means infinity)
%  minpts - minimum number of data points to perform fit
%      (set to 5 if smaller number if given)
%  te - EFW time in seconds (isGetDataLite time)
%  data - EFW data from pair in mV/m, should correspond to te
%  tp - Ephemeris time in seconds (isGetDataLite time)
%  ph - Ephemeris phase in degr (sun angle for s/c Y axis), should 
%      correspond to tp 
%  method - 0: c_efw_onesfit, Matlab routine by AIE
%           1: c_efw_spinfit_mx (default), BHN Fortran source obtained from KTH
%
% Output:
%  spinfit = [ts,ex,ey,offset,sdev0,sdev,iter,nout]
%  ts - time vector in seconds 
%  ex - E-field x-component in DSI coordinates (almost GSE)
%  ey - E-field y-component in DSI coordinates (almost GSE)
%  offset - mean value of input data
%  sdev0 - standard deviation in first fit. Has no meaning if METHOD=1
%  sdev - standard deviation in final fit
%  iter - number of iterations (one if OK at once)
%  nout - number of outliers removed
%
% This function chops up the time series in four second
% intervals, each of which are analysed by c_efw_onesfit.
% Spins with less then 90% of data are disregarded.
%
% Example: (assuming valid db exist, see isGetDataLite)
% To plot spinfits of data from probe pair 34 (obtained during
%    the first 10 minutes of the Donald Duck show on Swedish
%    national TV on Christmas Eve 2001), defining outliers as
%    points outside 3 standard deviations, allowing maximum 10
%    iterations for each spin, and demanding 20 data points for
%    the fit to be valid, do as follows:
%  [t34,e34] = isGetDataLite(db, [2002 12 24 14 00 00], 600, ...
%  'Cluster', '4', 'efw','E','p34','10Hz','hx');  
%  [tpha,pha] = isGetDataLite(db, [2002 12 24 14 00 00], 600, ...
%  'Cluster', '4', 'ephemeris','phase','','','');  
%  sp34 = c_efw_sfit(34,3,10,20,t34,e34,tpha,pha);
%  t0 = toepoch([2002 12 24 14 00 00]);
%  t = sp34(:,1) - t0;
%  ex = sp34(:,2);
%  ey = sp34(:,3);
%  plot(t,ex,'k',t,ey,'r');
%  xlabel('Time from 2002-12-24 14:00:00 [s]');
%  ylabel('E [mV/m]');
%  title('Cluster SC4 EFW Spin Fits (Ex black, Ey red)')
%
% See also C_EFW_ONESFIT, isGetDataLite
%
% $Id$

%
% Original version by Anders.Eriksson@irfu.se, 13 December 2002

error(nargchk(8,9,nargin))

if pair~=12 && pair~=32 && pair~=34, error('PAIR must be one of: 12, 32, 34'), end

% Set default method to BHN
if nargin < 9, method = 1; end

if method==1
	if exist('c_efw_spinfit_mx','file')~=3
		method = 0;
		disp('cannot find mex file, defaulting to Matlab code.')
	end
end	

% Turn off warnings for badly conditioned polynomial:
warning off;

% Chop up time interval
% We always start at 0,4,8.. secs, so that we have 
% the same timelines an all SC at 2,6,10... sec
tstart = fix(min(te)/4)*4;
tend = max(te);
n = floor((tend - tstart)/4) + 1;
spinfit = zeros(n,8);

% Guess the sampling frequency
sf = c_efw_fsample(te);

% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .98; 

n_gap = 0;

if method ==1
	fnterms = 3;
	te = torow(te);
	data = torow(data);
	x = zeros(n,9);
	phi = zeros(n,1);
	spinfit(:,[5 8]) = -1;

	tpha = tocolumn(tp);
	pha = tocolumn(ph);
	% Calcluate phase (in rad) at EFW sample times:
	pha = unwrap(pi*pha/180);
	% Find phase of given pair:
	if pair == 12, pha = pha + 3*pi/4;
	elseif pair == 32, pha = pha + pi/2;
	elseif pair == 34, pha = pha + pi/4;
    else error('probe pair must be one of 12, 32 or 34')
	end
end

% Do it:
for i=1:n
	t0 = tstart + (i-1)*4;
	eind = find((te >= t0) & (te < t0+4));
	pind = find((tp >= t0) & (tp < t0+4));
	
	% Clear NaNs
	eind(isnan(data(eind))) = [];
	
	% Check for data gaps inside one spin.
	if sf>0 && length(eind)<N_EMPTY*4*sf, eind = []; end
	  
	% Need to check if we have any data to fit.
	if ~isempty(eind) && ~isempty(pind)
		if method==1
			% Use Fortran version of spin fit
			[bad,x(i-n_gap,:),spinfit(i-n_gap,6),spinfit(i-n_gap,7),lim] = ...
				c_efw_spinfit_mx(fnterms,maxit,2*pi/4.0,te(eind)-te(eind(1)),data(eind));
			tsfit = t0 + 2;	
			pol = polyfit(tpha(pind),pha(pind),1);
			phi(i-n_gap) = polyval(pol,tsfit);
			spinfit(i-n_gap,8) = length(find(bad==1));
			spinfit(i-n_gap,1) = tsfit;
		else
			% Use Matlab version by AIE
			spinfit(i - n_gap,:) = c_efw_onesfit(pair,fout,maxit,minpts,te(eind),...
				data(eind),tp(pind),ph(pind));
		end
    else n_gap = n_gap + 1;
	end 
end  
spinfit = spinfit(1:n - n_gap, :);

if method==1 && ~isempty(spinfit)
        x = x(1:n - n_gap, :);
        phi = phi(1:n - n_gap);

		theta = atan2(x(:,3),x(:,2));
		rho = sqrt(x(:,2).^2 + x(:,3).^2);
		
        % Correct phase
		% - Because s/c is spinning upside down:
        spinfit(:,2) = -rho.*cos(phi + theta);
        spinfit(:,3) = rho.*sin(phi + theta);
        spinfit(:,4) = x(:,1);
end

irf_log('proc',sprintf('%d spins processed, %d gaps found',n,n_gap))
