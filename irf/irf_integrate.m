function xint=irf_integrate(x,tref,time_step)
%IRF_INTEGRATE  Integrate time series
%
% xint=irf_integrate(x,tref,time_step)
%   integrate time series. time steps that are larger
%   than 3times the time step are assumed to be data gaps.
%
%   x - time series  to integrate
%   tref - optional, integration starting time (optional) 
%        isdat epoch or [yyyy mm dd hh mm ss.ss]
%   time_step - optional, all time_steps larger than 3*time_step are
%   assumed data gaps, default is that time_step is the smallest value of
%   all time_steps of the time series

isinpTS = isa(x,'TSeries');
if isinpTS
  unitsTmp = x.units;
  siConvTmp = x.siConversion;
  epochTmp = EpochTT(x.time).epoch;
  epoch0 = epochTmp(1);
  ttemp = double(epochTmp-epoch0)*1e-9;
  datatemp = double(x.data);
  x = [ttemp, double(datatemp)];
end

dt=[0 ; diff(x(:,1))];
if nargin < 3 % estimate time step
    time_steps=diff(x(:,1));
    [~,ind_min]=min(time_steps);
    time_steps(ind_min)=[]; % remove the smallest time step in case some problems
    time_step=min(time_steps);
end
dt(dt>3*time_step)=0;
xint=x;
for j=2:size(xint,2)
  j_ok=find(~isnan(xint(:,j)));
  xint(j_ok,j)=cumsum(x(j_ok,j).*dt(j_ok),1);
end

if nargin>=2 % other origo for integration 
    if isa(tref,'GenericTimeArray')
      tt = tref.epochUnix;
      if isinpTS, tt = tt - EpochTT(epoch0).epochUnix; end
    elseif size(tref)==6
        tt=toepoch(tref);
    elseif size(tref)==1
        tt=tref;
    else
        errS = 'do not know how to treat TREF input';
        irf.log('critical',errS), error(errS)
    end
    if tt < x(1,1),tt=x(1,1);end % if tref before first data point,set it to time of first data point
    if tt > x(end,1),tt=x(end,1);end % if tref after laast data point,set it to time of last data point
    xint_ref=irf_resamp(xint,tt,'nearest');
    xint=irf_add(1,xint,-1,xint_ref);
end

if isinpTS
  xintd = xint(:, 2:end);
  xintt = EpochTT(epochTmp);
  xint = TSeries(xintt,xintd);
  xint.units = sprintf('(%s)*s',unitsTmp);
  xint.siConversion = [siConvTmp ' s'];
end