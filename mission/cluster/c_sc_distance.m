function [sc_distance]=c_sc_distance(tint,dt)
%C_SC_DISTANCE  estimate and plot the distance between spacecraft
%
%   [sc_distance]=c_sc_distance(tint,dt);
%
% INPUT
%       tint - time interval for which to plot distance
%               [tstart tend] in isdat epoch
%         dt - time step, if not given then one point for each orbit
%
% OUTPUT
%       sc_distance=[time dr1-dr2 dr1-dr3 dr1-dr4 dr2-dr3 dr2-dr4 dr3-dr4]
%
%   c_sc_distance; % plot for the whole cluster interval
%

narginchk(0,2)
if nargin==0
  help c_sc_distance;
  disp('******************************************************')
  disp('************ Running for the whole life time  ********')
  disp('******************************************************')
  tint=[toepoch([2000 01 01 0 0 0]) date2epoch(datenum(date))];
  dt=0;
elseif nargin==1 % definte tint as one orbit
  disp('************ One data point each orbit  ********')
  dt=0;
end

sc_distance=[];
if dt==0 % find the apogee times
  c_eval('data=getData(ClusterDB(''db:10''),tint(1),tint(2)-tint(1),?,''r'',''nosave'');R?=data{2};');
  c_eval('r?=irf_resamp(R?,R1);');
  c_eval('r!r?=irf_add(1,r!,-1,r?);',1:4,1:4);
  c_eval('r!r?=irf_abs(r!r?);',1:4,1:4);
  c_eval('r?=irf_abs(r?);');
  rsum=r1r1(:,1:2);
  c_eval('rsum(:,2)=rsum(:,2)+r!r?(:,5);',1:4,1:4);
  index_min_separation=[];
  for jj=2:(size(rsum,1)-1)
    if rsum(jj,2)<rsum(jj-1,2) && rsum(jj,2)<rsum(jj+1,2)
      index_min_separation=[index_min_separation jj];
    end
  end
  imin=index_min_separation;
  sc_distance=[r1r2(imin,[1 5]) r1r3(imin,5) r1r4(imin,5) r2r3(imin,5) r2r4(imin,5) r3r4(imin,5)];
else
  for tt=tint(1):dt:tint(2)
    c_eval('data=getData(ClusterDB,tt,dt,?,''r'',''nosave'');R?=data{2};');
    c_eval('r?=irf_resamp(R?,R1);');
    c_eval('r!r?=irf_add(1,r!,-1,r?);',1:4,1:4);
    c_eval('r!r?=irf_abs(r!r?);',1:4,1:4);
    rsum=r1r1(:,1:2);
    c_eval('rsum(:,2)=rsum(:,2)+r!r?(:,5);',1:4,1:4);
    [rsummin,imin]=min(rsum(:,2));
    sc_distance(end+1,:)=[r1r2(imin,[1 5]) r1r3(imin,5) r1r4(imin,5) r2r3(imin,5) r2r4(imin,5) r3r4(imin,5)];
  end
end

figure;
ccol=['k','r','g','b'];
plot(sc_distance(:,1),sc_distance(:,2),'-ko','MarkerFaceColor','r','MarkerSize',2);
title('Sc separation');
ylabel('dR [km]');
irf_pl_info(['c\_sc\_distance() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); % add information to the plot

grid on
hold on
isc=2;
for ic1=1:3
  for ic2=ic1+1:4
    c_eval('plot(sc_distance(:,1),sc_distance(:,isc),[''--'',ccol(?)],''LineWidth'',3);',ic2);
    c_eval('plot(sc_distance(:,1),sc_distance(:,isc),[''-o'',ccol(?)],''MarkerFaceColor'',ccol(?),''MarkerSize'',5);',ic1);
    isc=isc+1;
  end
end
irf_timeaxis(gca,'date');


