function flag=c_ri_get_event_data(time_interval,path_Events,path_Out, data_list, dt_interval)
%function flag=c_ri_get_event_data(time_interval,path_Events,path_Out, data_list, dt_interval)
%
%Input:
%   time_interval - [start_time end_time] in isdat_epoch
%   path_Events - where to look for Events files (files start with "E_"), alternative vector with times of efents in isdat epoch [t1 t2 t3 ..]
%   path_Out    - where to write the data, ends with '/'
%   data_list   - structure with values of which data to load, e.g. {'EFW_E','EFW_P'}
%   dt_interval - download +- dt_interval seconds around event (default 5s)
%
%Descrition of the function:
%   Loads the specified data for each event in the given time interval
%
%Description of variables:
% eventsw   - 1 event per row, 3 columns [start_time end_time s/c_mode(1-burst, 0-normal)]
%
%Written by Andris Vaivads Sep 2003 -

global AV_DEBUG; if isempty(AV_DEBUG), debug=0;else, debug=AV_DEBUG; end

%--------------------- the beginning --------------------------
sc_list=1:4; % get data for all 4 s/c
isdat_database='disco:10'; %
db = Mat_DbOpen(isdat_database);

if nargin==0, help c_ri_get_event_data; return; end
if nargin==1, path_Events=[pwd filesep];path_Out=[pwd  filesep];data_list={'EFW_P'};dt_interval=5;end
if nargin==2, path_Out=[pwd  filesep];data_list={'EFW_P'};dt_interval=5;end
if nargin==3, data_list={'EFW_P'};dt_interval=5;end
if nargin==4, dt_interval=5;end
default_cases={'EPH','FGM'};
data_list=[default_cases data_list ]; % ephemeris should be first


% construct time intervals to download, time intervals are in whole seconds
if isstr(path_Events)
  start_time=time_interval(1);
  end_time=time_interval(2);
  event_time_intervals=[]; % three columns [start_time end_time mode]; mode=0(normal), 1(burst)
  next_event_row=1;
  dir_list=dir([path_Events 'E_' '*.mat']);
  A_list=dir([path_Events 'A/Ap_*.mat']);
  for i_Event_file=1:size(dir_list,1)
    event_file = dir_list(i_Event_file).name;
    load([path_Events event_file]); % load time_of_events variable
    flag_time_within_interval=sign((end_time-time_of_events(:,1)).*(time_of_events(:,1)-start_time));
    ind_events=find(flag_time_within_interval == 1);
    if ind_events
      found_events=[floor(time_of_events(ind_events,1)-dt_interval) ceil(time_of_events(ind_events,1)+dt_interval) time_of_events(ind_events,5)];
      event_time_intervals(next_event_row+[0:size(found_events,1)-1],:)=found_events;
    end
  end
elseif isnumeric(path_Events)
  event_time_intervals=[floor(path_Events(:)-dt) ceil(path_Events(:)+dt)];
else
  error(' events not defined properly');
end


disp(['Found ' num2str(size(event_time_intervals,1)) ' events.']);

% clean up overlapping time intervals
events=event_time_intervals(1,:);i_final_events=1;
for i_event=2:size(event_time_intervals,1)
  if event_time_intervals(i_event,1)<=events(i_final_events,2)
    events(i_final_events,2)=event_time_intervals(i_event,2);
  else
    i_final_events=i_final_events+1;
    events(i_final_events,:)=event_time_intervals(i_event,:);
  end
end
disp(['From ' num2str(size(event_time_intervals,1)) ' events constructed ' num2str(size(events,1)) ' nonoverlapping events.']);

if exist('mWork.mat','file'), save -append mWork events, else, save mWork events; end

for i_event=1:size(events,1)
  start_time_epoch=events(i_event,1);start_time=fromepoch(start_time_epoch);
  end_time_epoch  =events(i_event,2);end_time=fromepoch(end_time_epoch);
  time_interval=[start_time_epoch end_time_epoch];
  Dt        =end_time_epoch-start_time_epoch;
  if debug, disp([num2str(i_event) '.event, ' R_datestring(start_time) ', dt=' num2str(Dt) 's.']);end
  % sc_mode estimate fast solution
  sc_mode=[];
  for i_a=1:size(A_list,1)
    a_file=A_list(i_a).name;
    if c_ri_timestr_within_intervall(a_file,start_time,end_time) == 1
      sc_mode=a_file(length(a_file)-4);
    end
  end
  if isempty(sc_mode), disp('do not know which mode satellites are running, assuming normal!');sc_mode='n';end
  if debug, disp(['sc_mode=' sc_mode]);end

  for i_data=1:length(data_list)
    switch data_list{i_data}
      case 'EPH' % get ephemeris R,V,A,ILAT,MLT, + (not implemented but necessary) satellite axis orientation
        file_prefix='F';
        file_name=[path_Out file_prefix deblank(R_datestring(start_time)) '_T' deblank(R_datestring(end_time)) '.mat'];
        for ic=sc_list
          if debug, disp(['Loading ephemeris s/c' num2str(ic)]);end
          [tlt,lt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'lt', ' ', ' ', ' ');
          [tmlt,mlt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'mlt', ' ', ' ', ' ');
          [tL,Lshell] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'l_shell', ' ', ' ', ' ');
          [tilat,ilat] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'inv_lat', ' ', ' ', ' ');if debug,disp('LT,MLT,L shell,ILAT...ready!');end
          [tr,r] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
          [tv,v] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
          [tA,A] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');if debug,disp('Position R, velocity V, Phase A ...ready!');end
          eval(av_ssub('A?=[double(tA) double(A)];',ic));
          eval(av_ssub('LT?=[double(tlt) double(lt)];MLT?=[double(tmlt) double(mlt)];L?=[double(tL) double(Lshell)];ILAT?=[double(tilat) double(ilat)];R?=[double(tr) double(r)''];V?=[double(tv) double(v)''];',ic));clear tlt tmlt tL tilat lt mlt Lshell ilat tr r tv v;
          if exist(file_name,'file'), flag_append='-append';else, flag_append='';end
          stric=num2str(ic);
          save(file_name,['A' stric],['L' stric],['LT' stric],['MLT' stric],['ILAT' stric],['R' stric],['V' stric],flag_append);
          if debug, disp(['saving ephemeris for sc' num2str(ic)]);end
        end

      case 'FGM'
        file_prefix='F';
        file_name=[path_Out file_prefix deblank(R_datestring(start_time)) '_T' deblank(R_datestring(end_time)) '.mat'];
        [B1,B2,B3,B4]=c_get_bfgm(time_interval);
        for ic=sc_list,eval(av_ssub('dB?=c_gse2dsc(B?,?);',ic)),end
        if exist(file_name,'file'), flag_append='-append';else, flag_append='';end
        save(file_name,'B1','B2','B3','B4','dB1','dB2','dB3','dB4',flag_append);
        if debug, disp('saving B1 B2 B3 B4');end

      case 'EFW_P'
        file_prefix='F';
        file_name=[path_Out file_prefix deblank(R_datestring(start_time)) '_T' deblank(R_datestring(end_time)) '.mat'];
        EFW_P=c_isdat_get_EFW(time_interval,[],[],sc_mode,1:4,db,'P');
        P1=EFW_P{1};P2=EFW_P{2};P3=EFW_P{3};P4=EFW_P{4};
        if exist(file_name,'file'), flag_append='-append';else, flag_append='';end
        save(file_name,'P1','P2','P3','P4',flag_append);        if debug, disp(['saving ' flag_append ' P1,P2,P3,P4 ->' file_name]);end

      case 'EFW_E'
        file_prefix='F';
        file_name=[path_Out file_prefix deblank(R_datestring(start_time)) '_T' deblank(R_datestring(end_time)) '.mat'];
        deg=20; % the minimum elevation of B with respect to the spin plane when E.B=0 is used for spin axis E
        for ic=sc_list
          eval(av_ssub('wE?=c_isdat_get_EFW(time_interval,[],[],sc_mode,?,db,''wE'');',ic));
          eval(av_ssub('dE?=c_despin(wE?,?,''efw'');',ic)),
          eval(av_ssub('deg=20;[dE?,d?]=av_ed(dE?,dB?,deg);E?=c_gse2dsc(dE?,[dE?(1,1) ?],-1);indzero=find(abs(d?)<deg);E?(indzero,4)=0;',ic));
          eval(av_ssub('ExB?=av_e_vxb(E?,B?,-1);',ic));
          if exist(file_name,'file'), flag_append='-append';else, flag_append='';end
          eval(av_ssub('save(file_name,''wE?'',''dE?'',''d?'',''E?'',''ExB?'',flag_append);',ic));
          if debug, disp(['saving wE dE d E ExB for sc' num2str(ic)]);end
        end
    end
  end
end
Mat_DbClose(db);

flag=0; % output not yet implemented

