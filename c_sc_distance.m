function [d1]=c_sc_distance(tint,dt);
%C_SC_DISTANCE  estimate and plot the distance between spacecraft
%
%   [d1]=c_sc_distance(tint,dt);
%
% INPUT 
%       tint - time interval for which to plot distance 
%               [tstart tend] in isdat epoch
%         dt - time step, if not given then one point for each orbit 
% 
%   c_sc_distance; % plot for the whole cluster interval 
%
% $Id$

error(nargchk(0,2,nargin))
if nargin==0,
    tint=[toepoch([2000 01 01 0 0 0]) date2epoch(datenum(date))];
elseif nargin==1, % definte tint as one orbit
    tint=52*3600;
end

sc_distance=[];
for tt=tint(1):dt:tint(2),
    c_get_batch(tt,dt,'sp','./','vars',{'r'});
    load mR;
    c_eval('r?=irf_resamp(r?,r1);');
    c_eval('r!r?=irf_add(1,r!,-1,r?);',1:4,1:4));
    c_eval('r!r?=irf_abs(r!r?);',1:4,1:4));
    rsum=r1r1(:,1:2);
    c_eval('rsum(:,2)=rsum(:,2)+r!r?(:,5);',1:4,1:4);
    [rsummin,imin]=min(rsum(:,2));
    sc_distance(end+1,:)=[r1r2(imin,[1 5]) r1r3(imin,5) r1r4(imin,5) r2r3(imin,5) r2r4(imin,5) r3r4(imin,5)];
end

ccol=['k','r','g','b'];
h=irf_plot(sc_disctance(:,1:2),'-ko','MarkerFaceColor','r','MarkerSize',2);
title('Sc separation');
ylabel('dR [km]');
irf_pl_info(['c\_sc\_distance() ' datestr(now)]); % add information to the plot

grid on
hold on
isc=3;
for ic1=1:4,
 for ic2=ic1+1:4,
  c_eval('plot(sc_distance(:,1),sc_distance(:,isc3),[''-'',ccol(?),''o''],''MarkerFaceColor'',ccol(!),''MarkerSize'',2);',ic1,ic2));
 end
end
add_timeaxis(gca,'date');


