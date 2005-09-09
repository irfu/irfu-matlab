function xint=av_integrate(x,tref)
%IRF_INTEGRATE  Integrate time series
%
% xint=irf_integrate(x,tref)
%   integrate time series. time steps that are larger
%   than 3times the smallest time step are assumed to be data gaps
%
%   x - time series  to integrate
%   tref - optional, integration starting time (optional) 
%        isdat epoch or [yyyy mm dd hh mm ss.ss]
%
% $Id$

dt=[0 ; diff(x(:,1))];
time_step=min(diff(x(:,1)));
data_gaps=find(diff(x(:,1))>3*time_step);
dt(data_gaps)=0;
xint=x;
xint(:,2:end)=cumsum(x(:,2:end).*repmat(dt,1,size(x,2)-1),1);

if nargin==2, % other origo for integration 
    if size(tref==6),
        tt=toepoch(tref);
    elseif size(tref==1);
        tt=tref;
    else
        irf_log('proc','do not know how to treat <tref> input')
    end
    if tt < x(1,1),tt=x(1,1);end % if tref before first data point,set it to time ofo first data point
    if tt > x(end,1),tt=x(end,1);end % if tref after laast data point,set it to time ofo last data point
    xint_ref=irf_resamp(xint,tt);
    xint=irf_add(1,xint,-1,xint_ref);
end

