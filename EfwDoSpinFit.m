function spinfit = EfwDoSpinFit(pair,fout,maxit,minpts,te,data,tp,ph,method)
% EfwDoSpinFit -- will produce spin fit data values (ex, ey)
%   of EFW data frome given probe pair at 4 second resolution.
%
% Input:
%  pair - probe pair used (12, 34)
%  fout - minimum fraction of fit standard deviation that defines an outlier 
%         (zero means no removal of outliers)
%  maxit - maximum number of iterations (zero means infinity)
%  minpts - minimum number of data points to perform fit
%      (set to 5 if smaller number if given)
%  te - EFW time in seconds (isGetDataLite time)
%  data - EFW data from pair in mV/m, should correspond to te
%  tp - Ephemeris time in seconds (isGetDataLite time)
%  ph - Ephemeris phase in degr (sun angle for s/c Y axis), should 
%      correspond to tp 
%  method - 0: EfwDoOneSpinFit (default, Matlab routine by AIE), 
%           1: spinfit_mx MEX file, Fortran source obtained from PAL, KTH.
%
% Output:
%  spinfit = [ts,ex,ey,offset,sdev0,sdev,iter,nout]
%  ts - time vector in seconds 
%  ex - E-field x-component in DSI coordinates (almost GSE)
%  ey - E-field y-component in DSI coordinates (almost GSE)
%  offset - mean value of input data
%  sdev0 - standard deviation in first fit
%  sdev - standard deviation in final fit
%  iter - number of iterations (one if OK at once)
%  nout - number of outliers removed
%
% This function chops up the time series in four second
% intervals, each of which are analysed by EfwDoOneSpinFit.
% Spins with less then 75% of data are disregarded.
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
%  sp34 = EfwDoSpinFit(34,3,10,20,t34,e34,tpha,pha);
%  t0 = toepoch([2002 12 24 14 00 00]);
%  t = sp34(:,1) - t0;
%  ex = sp34(:,2);
%  ey = sp34(:,3);
%  plot(t,ex,'k',t,ey,'r');
%  xlabel('Time from 2002-12-24 14:00:00 [s]');
%  ylabel('E [mV/m]');
%  title('Cluster SC4 EFW Spin Fits (Ex black, Ey red)')
%
% See also: EfwDoOneSpinFit, isGetDataLite
%
% $Id$
%
% Original version by Anders.Eriksson@irfu.se, 13 December 2002

if nargin < 9, method = 0; end
if method==1
	if exist('spinfit_mx')~=3
		method = 0;
		disp('cannot find mex file, defaulting to Matlab code.')
	end
end	

% Turn off warnings for badly conditioned polynomial:
warning off;

% Chop up time interval:
tstart = min(te);
tend = max(te);
n = floor((tend - tstart)/4);
spinfit = zeros(n,8);

% guess the sampling frequency
sf = length(data)/(tend - tstart);

if sf<1.3*25 & sf>.7*25, sf = 25;
elseif sf<1.3*450 & sf>.7*450, sf = 450;
else
	disp('cannot guess sampling frequency')
	sf = 0;
end

% N_EMPTY .75 means that we use only sping with more then 75% points.
N_EMPTY = .75; 

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
	elseif pair == 34, pha = pha + pi/4;
	else, error('probe pair must be one of 12 or 34')
	end
end

% Do it:
for i=1:n
	t0 = tstart + (i-1)*4;
	eind = find((te >= t0) & (te < t0+4));
	pind = find((tp >= t0) & (tp < t0+4));

	% check for data gaps inside one spin.
	if sf>0 & length(eind)<N_EMPTY*4*sf, eind = []; end
	  
	% wee need to check if we have any data to fit.
	if ~isempty(eind) & ~isempty(pind)
		if method==1
			%we use Fortran version of spin fit
			[bad,x(i-n_gap,:),spinfit(i-n_gap,6),spinfit(i-n_gap,7),lim] = ...
				spinfit_mx(fnterms,maxit,2*pi/4.0,te(eind),data(eind));
			
			tsfit = mean(te(eind));	
			pol = polyfit(tpha(pind),pha(pind),1);
			phi(i-n_gap) = polyval(pol,tsfit);

			spinfit(i-n_gap,1) = tsfit;
		else
			%we use Matlab version by AIE
			spinfit(i - n_gap,:) = EfwDoOneSpinFit(pair,fout,maxit,minpts,te(eind), ...
				data(eind),tp(pind),ph(pind));
		end
	else, n_gap = n_gap + 1;
	end 
end  
spinfit = spinfit(1:n - n_gap, :);

if method==1
        x = x(1:n - n_gap, :);
        phi = phi(1:n - n_gap);

		theta = atan2(x(:,3),x(:,2));
		rho = sqrt(x(:,2).^2 + x(:,3).^2);
		
        %correct phase
		% - Because s/c is spinning upside down:
        spinfit(:,2) = -rho.*cos(phi + theta);
        spinfit(:,3) = rho.*sin(phi + theta);
        spinfit(:,4) = x(:,1);
end

c_log('proc',sprintf('%d spins processed, %d gaps found',n,n_gap))
