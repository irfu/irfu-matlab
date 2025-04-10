function [soloPosition, bepiPosition, mercuryPosition, venusPosition, marsPosition, earthPosition]= plot_orbit(tint,ttime,varargin);
% A function to plot and get the position of Solar Orbiter [and Mars, Venus, Earth, Mercury, BepiColombo] for a given time interval
% and a given time. Output: TSeries with position(s).
%
% Example
% tint = irf.tint('2020-02-11T00:00:00Z/2020-12-01T15:00:00Z'); % time interval in TT2000 UTC
% Ttime=[2020 07 12]; %Specific time of interest
% [Solo,Be, Me,V,Ma,E]=solo.plot_orbit(tint,Ttime,['Bepi'],['Mercury'],['Venus'],['Earth'],['Mars']);

frame='ECLIPJ2000';
t_step = 60*60; %3600 sec

% Find the most recent SolO metakernel file and load it
sharedPath = '/share/SPICE/'; % SPICE kernels for different missions are found in this folder on IRFU servers
dirsBep = dir([sharedPath,'BEPICOLOMBO/kernels/mk/bc_plan_v*_*.tm']);
if size(dirsBep, 1) > 1
  % Multiple kernels could be found if executing this script at the same time as syncing new kernel files
  error('Found multiple metakernels, please check your folder.');
end
kernelFileBep = [dirsBep.folder, filesep, dirsBep.name];
cspice_furnsh(kernelFileBep);

% Find the most recent BepiColombo metakernel file and load it %Ihope it is
% ok to load two metakernels...
sharedPath = '/share/SPICE/'; % SPICE kernels for different missions are found in this folder on IRFU servers
dirs = dir([sharedPath,'Solar-Orbiter/kernels/mk/*pred-mk_v*.tm']);
if size(dirs, 1) > 1
  % Multiple kernels could be found if executing this script at the same time as syncing new kernel files
  error('Found multiple metakernels, please check your folder.');
end
kernelFile = [dirs.folder, filesep, dirs.name];
cspice_furnsh(kernelFile);





% Compute et (SPICE ephemeries time, make use of input in TT2000)
et = tint.start.tts:t_step:tint.stop.tts;

% The position of Solar Orbiter in units of km
pos_solo = cspice_spkpos('solo', et, frame, 'none', 'Sun');

% Convert to utc and then to TT2000 in TSeries object
utc_tmp = cspice_et2utc(et, 'ISOC', 0);

% Note pos' since it is returned as 3xN but TSeries expects Nx3 (where N is number of records).
soloPosition = irf.ts_vec_xyz(EpochTT(utc_tmp), pos_solo'/1.496e8);
soloPosition.units = 'km'; % Add some metadata information (read when plotting in irf_plot)
soloPosition.name='Solar Orbiter';
soloPosition.coordinateSystem=frame;

f1=figure(1);clf
f1.Position=[10 1800 1100 450];

%Plot the positio of Solar Orbiter
irf_subplot(1,2,1)
plot(soloPosition.data(:,1),soloPosition.data(:,2),'k-'); hold on; grid on
ylabel('Y [AU]'); xlabel('X [AU]')

%Plot the Positions at a single time (Ttime)
t=irf_time(soloPosition.time,'epochtt>date');
Ttime=irf_time(ttime,'vector>date');
mini=find(min(abs(t-Ttime))==abs(t-Ttime));
plot(soloPosition.data(mini,1),soloPosition.data(mini,2),'ks');
timetxt=datestr(irf_time(ttime,'vector>date'),'YYYY-mm-dd');
text(soloPosition.data(mini,1)*1.1,soloPosition.data(mini,2)*1.1,['SolO',' ',timetxt], 'fontsize', 8, 'Color','k')

%Plot the postion(s) of various planets
if any(ismember(varargin,'Bepi'))
  pos_bepi = cspice_spkpos('BEPICOLOMBO MMO', et, frame, 'none', 'Sun');
  bepiPosition= irf.ts_vec_xyz(EpochTT(utc_tmp), pos_bepi'/1.496e8);
  bepiPosition.units = 'AU'; % Add some metadata information (read when plotting in irf_plot)
  bepiPosition.coordinateSystem=frame;
  bepiPosition.name='BepiColombo MMO';
  plot(bepiPosition.data(:,1),bepiPosition.data(:,2),'-','Color',[0.7 0.7 0])
  plot(bepiPosition.data(mini,1),bepiPosition.data(mini,2),'s','Color',[0.4 0.4 0])
  text(bepiPosition.data(mini,1)*1.1,bepiPosition.data(mini,2)*1.1,'BepiColombo', 'fontsize', 8, 'Color',[0.4 0.4 0])
end

if any(ismember(varargin,'Mercury'))
  pos_mercury = cspice_spkpos('mercury', et, frame, 'none', 'Sun');
  mercuryPosition= irf.ts_vec_xyz(EpochTT(utc_tmp), pos_mercury'/1.496e8);
  mercuryPosition.units = 'AU'; % Add some metadata information (read when plotting in irf_plot)
  mercuryPosition.coordinateSystem=frame;
  mercuryPosition.name='Mercury';
  plot(mercuryPosition.data(:,1),mercuryPosition.data(:,2),'-','Color',[0.4 0.4 0])
  plot(mercuryPosition.data(mini,1),mercuryPosition.data(mini,2),'s','Color',[0.4 0.4 0])
  text(mercuryPosition.data(mini,1)*1.1,mercuryPosition.data(mini,2)*1.1,'Mercury', 'fontsize', 8, 'Color',[0.4 0.4 0])
end

if any(ismember(varargin,'Venus'))
  pos_venus = cspice_spkpos('venus', et, frame, 'none', 'Sun');
  venusPosition= irf.ts_vec_xyz(EpochTT(utc_tmp), pos_venus'/1.496e8);
  venusPosition.units = 'AU'; % Add some metadata information (read when plotting in irf_plot)
  venusPosition.coordinateSystem=frame;
  venusPosition.name='Venus';
  plot(venusPosition.data(:,1),venusPosition.data(:,2),'-','Color',[0.8 0.6 0])
  plot(venusPosition.data(mini,1),venusPosition.data(mini,2),'s','Color',[0.8 0.6 0])
  text(venusPosition.data(mini,1)*1.1,venusPosition.data(mini,2)*1.1,'Venus', 'fontsize', 8, 'Color',[0.8 0.6 0])
end

if any(ismember(varargin,'Mars'))
  pos_mars = cspice_spkpos('mars', et, frame, 'none', 'Sun');
  marsPosition= irf.ts_vec_xyz(EpochTT(utc_tmp), pos_mars'/1.496e8);
  marsPosition.units = 'AU'; % Add some metadata information (read when plotting in irf_plot)
  marsPosition.coordinateSystem=frame;
  marsPosition.name='Mars';
  plot(marsPosition.data(:,1),marsPosition.data(:,2),'r-')
  plot(marsPosition.data(mini,1),marsPosition.data(mini,2),'rs')
  text(marsPosition.data(mini,1)*1.1,marsPosition.data(mini,2)*1.1,'Mars', 'fontsize', 8, 'Color','r')
end

if any(ismember(varargin,'Earth'))
  pos_earth = cspice_spkpos('earth', et, frame, 'none', 'Sun');
  earthPosition= irf.ts_vec_xyz(EpochTT(utc_tmp), pos_earth'/1.496e8);
  earthPosition.units = 'AU'; % Add some metadata information (read when plotting in irf_plot)
  earthPosition.coordinateSystem=frame;
  earthPosition.name='Earth';
  plot(earthPosition.data(:,1),earthPosition.data(:,2),'b-')
  plot(earthPosition.data(mini,1),earthPosition.data(mini,2),'bs')
  text(earthPosition.data(mini,1)*1.1,earthPosition.data(mini,2)*1.1,'Earth', 'fontsize', 8, 'Color','b')
end

plot(0,0,'yo','MarkerFaceColor','y')
axis('equal')
axis('square')
ylim([-1.7,1.7])
xlim([-1.7,1.7])
text(0.03,0.03,frame, 'Units', 'Normalized', 'fontsize', 8)
str=[datestr(irf_time(tint(1),'epochtt>date'),'YYYY-mm-dd'),' - ',datestr(irf_time(tint(end),'epochtt>date'),'YYYY-mm-dd')];
title(str)

irf_subplot(1,2,2)
irf_plot(soloPosition.abs,'k-'); hold on
irf_plot(soloPosition(mini).abs,'ks')
ylabel('Heliocentric distance [AU]')
xlabel('UTC')
axis('square')
ylim([0,1.1])

%irf_print_fig(f1, 'soloposition','png')




