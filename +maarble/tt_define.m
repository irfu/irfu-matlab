% Find time intervals when Cluster and THEMIS are satisfying MAARBLE 
% requirements. The general datasets will include all the data for 
% 3 RE< R <10 RE 
% magnetic latitudes limited to 60 deg (to avoid the auroral zone).

% This software was developed as part of the MAARBLE (Monitoring,
% Analyzing and Assessing Radiation Belt Energization and Loss)
% collaborative research project which has received funding from the
% European Community's Seventh Framework Programme (FP7-SPACE-2011-1)
% under grant agreement n. 284520.

%% Cluster
Units=irf_units; %#ok<NASGU>
tintIso='2001-01-01T00:00:00.000Z/2013-01-01T00:00:00.000Z';
tint=irf_time(tintIso,'utc>tint'); %#ok<NASGU>
disp(['Using time interval: ' tintIso]);

disp('Loading Cluster 1-min positions');
load /data/caalocal/mR_1min; % load Cluster positions 1min resolution
c_eval('R?SM = irf.geocentric_coordinate_transformation(R?,''GSE>SM'');')
c_eval('R?=R?SM;clear R?SM;')
c_eval('izero=find(R?(:,1)==0);R?(izero,:)=[];');
c_eval('R?=irf_tlim(R?,tint);');
c_eval('R?=irf_abs(R?);');
c_eval('RRE?=irf_tappl(R?(:,[1 5]),''*Units.km/Units.RE'');');

disp('Calculating Cluster mlat, saving to matMlat');
c_eval('mlat?=[R?(:,1) asin(R?(:,4)./R?(:,5))*180/pi];');
save /data/caalocal/matMlat mlat1 mlat2 mlat3 mlat4
c_eval('clear R?;');

%load /data/caa/CAA/matMlat
disp('Finding when Cluster satisfies MAARBLE conditions');
tStep=median(diff(RRE1(:,1))); %#ok<NASGU> % time step 
minR=3; %#ok<NASGU> % minimum distance from Earth 
maxR=10; %#ok<NASGU>  % minimum distance from Earth 
maxMlat=60; % maximum magnetic latitude

% maarble definition
ttLabel='MAARBLE';
ttTitle='Cluster ? inside MAARBLE area, 3RE<R<10RE,mlat<60deg';
c_eval('imaarble?=(RRE?(:,2)>minR & RRE?(:,2)<maxR & abs(mlat?(:,2)) < maxMlat);')
% define intervals
c_eval('indstart?=find(diff([0 imaarble?(:)'']) == 1);');
c_eval('indend?=find(diff([imaarble?(:)'' 0]) == -1);');
c_eval(['clear tt_C?_in_' ttLabel])
c_eval(['tt_C?_in_' ttLabel '=irf.TimeTable;'])
c_eval(['tt_C?_in_' ttLabel '.Header={''' ttTitle '''};']);
c_eval(['tt_C?_in_' ttLabel '.TimeInterval=[RRE?(indstart?,1)-tStep/2 RRE?(indend?,1)+tStep/2];'])
c_eval(['disp(''Created time table: tt_C?_in_' ttLabel ''');']);
c_eval(['tt_C?_in_' ttLabel '=remove(tt_C?_in_' ttLabel ',find(diff(tt_C?_in_' ttLabel '.TimeInterval,1,2)<10*60));']);
y=irf_ask('Shall I save the time tables to IRF yes/no? [%]','y','no');  
if strcmp(y,'yes')
	c_eval(['irf.tt(tt_C?_in_' ttLabel ',''write_IRF'',''C?_in_' ttLabel ''');'])
end

%% THEMIS
Units=irf_units;
tintIso='2007-01-01T00:00:00.000Z/2012-01-01T00:00:00.000Z';
tint=irf_time(tintIso,'utc>tint');
disp(['Using time interval: ' tintIso]);

disp('Loading THEMIS 1-min positions');
load /data/caalocal/THEMIS/mRth.mat; % load THEMIS positions 1min resolution

c_eval('izero=find(Rth?(:,1)==0);Rth?(izero,:)=[];','abcde');
c_eval('Rth?=irf_tlim(Rth?,tint);Rth?=irf_abs(Rth?);','abcde');
c_eval('RRE?=irf_tappl(Rth?(:,[1 5]),''*Units.km/Units.RE'');','abcde');

disp('Finding when THEMIS satisfies MAARBLE conditions');
tStep=median(diff(RREa(:,1))); % time step
minR=3;  % minimum distance from Earth 
maxR=10;  % minimum distance from Earth 

% maarble definition
ttLabel='MAARBLE';
ttTitle='THEMIS ? inside MAARBLE area, 3RE<R<10RE';
c_eval('imaarble?=(RRE?(:,2)>minR & RRE?(:,2)<maxR);','abcde')
% define intervals
c_eval('indstart?=find(diff([0 imaarble?(:)'']) == 1);','abcde');
c_eval('indend?=find(diff([imaarble?(:)'' 0]) == -1);','abcde');
c_eval(['clear tt_th?_in_' ttLabel],'abcde')
c_eval(['tt_th?_in_' ttLabel '=irf.TimeTable;'],'abcde','abcde')
c_eval(['tt_th?_in_' ttLabel '.Header={''' ttTitle '''};'],'abcde');
c_eval(['tt_th?_in_' ttLabel '.TimeInterval=[RRE?(indstart?,1)-tStep/2 RRE?(indend?,1)+tStep/2];'],'abcde')
c_eval(['disp(''Created time table: tt_th?_in_' ttLabel ''');'],'abcde');
c_eval(['tt_th?_in_' ttLabel '=remove(tt_th?_in_' ttLabel ',find(diff(tt_th?_in_' ttLabel '.TimeInterval,1,2)<10*60));'],'abcde');
y=irf_ask('Shall I save the time tables to IRF yes/no? [%]','y','no');  
if strcmp(y,'yes')
	c_eval(['irf.tt(tt_th?_in_' ttLabel ',''write_IRF'',''th?_in_' ttLabel ''');'],'abcde')
end
