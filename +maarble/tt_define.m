% find time intervals when Cluster is satisfying MAARBLE requirements
% The general datasets will include all the data for 
% 3 RE< R <10 RE 
% magnetic latitudes limited to 60 deg (to avoid the auroral zone).
%

Units=irf_units;
workDir=tempname;
mkdir(workDir);
cd(workDir);
tintIso='2001-01-01T00:00:00.000Z/2012-01-01T00:00:00.000Z';
tint=irf_time(tintiso,'iso2tint');
disp(['Using time interval: ' tintIso]);

disp('Loading Cluster 1-min positions');
load /data/caa/CAA/mR_1min; % load Cluster positions 1min resolution
c_eval('izero=find(R?(:,1)==0);R?(izero,:)=[];');
c_eval('R?=irf_tlim(R?,tint);');
c_eval('R?=irf_abs(R?);');
c_eval('RRE?=irf_tappl(R?(:,[1 5]),''*Units.km/Units.RE'');');
c_eval('clear R?;');

disp('Reading Cluster mlat');
c_eval('ILAT?=local.c_read(''Invar_Lat__C?_JP_PMP'',tint);');
c_eval('ILAT?(diff(ILAT?(:,1))==0,:)=[];')
c_eval('mlat?=irf_resamp(ILAT?,RRE?);');

disp('Finding when Cluster satisfies MAARBLE conditions');
tStep=median(diff(RRE1(:,1))); % time step
minR=3;  % minimum distance from Earth 
maxR=10;  % minimum distance from Earth 
maxMlt=60; % maximum magnetic latitude

% maarble definition
ttLabel='MAARBLE';
ttTitle='Cluster ? inside MAARBLE area, 3RE<R<10RE,mlat<60deg';
c_eval('imaarble?=(abs(RRE?(:,2))>minR & abs(RRE?(:,2))<maxR & mlat?(:,2) < maxMlt);')
% tailbox definition
% c_eval('itailbox?=(abs(RRE?(:,3))<tailBoxDY & abs(RRE?(:,4))<tailBoxDZ & RRE?(:,2)<tailBoxX);',sclist)
% define intervals for tailbox
c_eval('indstart?=find(diff([0 imaarble?(:)'']) == 1);');
c_eval('indend?=find(diff([imaarble?(:)'' 0]) == -1);');
c_eval(['tt_C?_in_' ttLabel '.start=RRE?(indstart?,1)-tStep/2;'])
c_eval(['tt_C?_in_' ttLabel '.end=RRE?(indend?,1)+tStep/2;'])
c_eval(['tt_C?_in_' ttLabel '.description=[''' ttTitle '''];']);
c_eval(['Created time table: ' 'tt_C?_in_' ttLabel ]);

y=irf_ask('Shall I save the time tables to IRF yes/no? [%]','y','no');  
if strcmp(y,'yes'),
	c_eval(['irf_tt(tt_C?_in_' ttLabel ',''write_IRF'',''C?_in_' ttLabel ''');'])
end