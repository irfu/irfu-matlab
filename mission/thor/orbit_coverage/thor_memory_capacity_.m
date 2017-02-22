%% Get omni data
tint = irf.tint('2004-01-01T00:00:00.00Z/2007-12-01T00:00:00.00Z');
%tint = irf.tint('2006-01-01T00:00:00.00Z/2009-12-01T00:00:00.00Z');
%tint = irf.tint('2004-01-01T00:00:00.00Z',60*60*20*30);
% omni_orig = thor_get_omni(tint,'Bx,By,Bz,bsnx,P,Ma,n,V,T'); % ~90 s/year
% save([irf('path') '/mission/thor/orbit_coverage/omni_data/omni2004-2008'],'omni_orig')
load([irf('path') '/mission/thor/orbit_coverage/omni_data/omni2004-2008.mat'])
omni = omni_orig;

%% Make TSeries of omni data and oversample NaN values
t_tmp = omni(:,1);
time = irf_time(t_tmp,'epoch>epochtt');

tsOMNI = irf.ts_scalar(time,omni(:,2:end));
isNaN_Bz = find(isnan(omni(:,4))); % Bz
isNaN_BSNX = find(isnan(omni(:,5))); % Bz
isNaN_Dp = find(isnan(omni(:,6))); % Dp
isNaN = unique([isNaN_Dp;isNaN_BSNX;isNaN_Bz]);
allInd = 1:time.length;
notNaN = setdiff(allInd,isNaN);
%negM = find((omni(:,7)<0));
%largeM = find((omni(:,7)>0));

tsB = irf.ts_vec_xyz(time(notNaN),omni(notNaN,2:4));
  tsB.units = 'nT';
  tsB.name = 'B';    
  tsB = tsB.resample(time);
tsBSNX = irf.ts_scalar(time(notNaN),omni(notNaN,5));
  tsBSNX.units = 'RE';
  tsBSNX.name = 'BSNX';  
  tsBSNX = tsBSNX.resample(time);
tsDp = irf.ts_scalar(time(notNaN),omni(notNaN,6));
  tsDp.units = 'nPa';
  tsDp.name = 'Dp';   
  tsDp = tsDp.resample(time);
tsM = irf.ts_scalar(time(notNaN),omni(notNaN,7));
  tsM.units = '';
  tsM.name = 'Ma';   
  tsM = tsM.resample(time);
tsN = irf.ts_scalar(time(notNaN),omni(notNaN,8));
  tsN.units = 'cc';
  tsN.name = 'n';   
  tsN = tsN.resample(time);
tsV = irf.ts_scalar(time(notNaN),omni(notNaN,9));
  tsV.units = 'km/s';
  tsV.name = 'V';     
  tsV = tsV.resample(time); 
tsT = irf.ts_scalar(time(notNaN),omni(notNaN,10));
  tsT.units = 'K?';
  tsT.name = 'T';    
  tsT = tsT.resample(time); 

%% Load THOR orbit
datastore('spice','dir','/Users/andris/calc/spice');
units = irf_units;

resampleKernels = 1; 
rTHOR_orig = thor_orbit('new1a.bsp',1*3600);
rTHOR = rTHOR_orig;
if resampleKernels
  newTime = rTHOR.time.start:60:rTHOR.time.stop; % 1 min intervals
  tmpR = rTHOR.resample(newTime);
  rTHOR = tmpR;
end

%% Resample THOR to OMNI timeline
tShift = time.start - rTHOR.time.start;
newTime = rTHOR.time + tShift; % shift the time of THOR to bsnx's time
rTHOR = irf.ts_vec_xyz(newTime,rTHOR.data); rTHOR.name = 'THOR orbit'; rTHOR.units = 'km';  

% cut the end of the TSeries
tsB    =    tsB.tlim(rTHOR.time);
tsBSNX = tsBSNX.tlim(rTHOR.time);
tsDp   =   tsDp.tlim(rTHOR.time);
tsM    =    tsM.tlim(rTHOR.time);
tsV    =    tsV.tlim(rTHOR.time);
tsN    =    tsN.tlim(rTHOR.time);
tsT    =    tsT.tlim(rTHOR.time);

rTHOR = rTHOR.resample(tsB); % upsample orbit times to OMNI timeline, 1 min

%% Define THOR phases, based on original timeline, add timeshift
tintTHOR = rTHOR_orig.time([1 end])+tShift;
endCalibrationPhase = irf_time('2026-12-31T00:00:12.994827148Z','utc>epochTT')+tShift; % 6 month commissioning phase
endPhase1 = irf_time('2027-12-10T12:00:12.994827148Z','utc>epochTT')+tShift;
endPhase2 = irf_time('2028-12-07T21:10:28.069581787Z','utc>epochTT')+tShift;
calibrationPhase = EpochTT([tintTHOR.start.utc; endCalibrationPhase.utc]);
%tintPhase1 = EpochTT([endCalibrationPhase.utc;endPhase1.utc]);
tintPhase1 = EpochTT([tintTHOR.start.utc;endPhase1.utc]);
tintPhase2 = EpochTT([endPhase1.utc; endPhase2.utc]);
tintPhase3 = EpochTT([endPhase2.utc; tintTHOR.stop.utc]);
tintSciencePhase = EpochTT([tintPhase1.start.utc; tintPhase3.stop.utc]);

%% Find local minima to divide complete trajectory into separate orbits
[valPerigee,isPerigee] = findpeaks(-rTHOR.abs.data*1e3/units.RE,'MinPeakProminence',3); %irf_plot({rTHOR.abs*1e3/units.RE},'comp'); hold on; irf_plot(rTHOR(isPerigee).abs*1e3/units.RE,'*')
tPerigee = rTHOR.time(isPerigee);
indOrbitsTPR = false(size(tPerigee));
indOrbitsTPR(rTHOR.data(isPerigee,1)<0 & ...
	abs(atan2d(rTHOR.data(isPerigee,2),-rTHOR.data(isPerigee,1)))<45) = true;

%irf_plot({rTHOR.abs*1e3/units.RE},'comp'); irf_plot(rTHOR(isPerigee).abs*1e3/units.RE,'*')

%% Check which KSR THOR is in
iKSR_orig = thor_in_ksrs(rTHOR,tsB,tsDp,tsM,tsBSNX);
iKSR = iKSR_orig; % iKSR can be manipulated later
FOM = iKSR;
FOM.data=zeros(size(iKSR.data));

%% Check quality factor of bowshock crossings
ind = 6; 
bsCrossing = find((iKSR.data)==ind); allInd = 1:iKSR.length; 
noCrossing = setdiff(allInd,bsCrossing);
Rout = 15;
QR = thor_QR(tsBSNX,Rout);
QV = thor_QV(tsBSNX,rTHOR);
[QBpar ,ShockNormalAngle] = thor_QB(tsBSNX,rTHOR,tsB);
[QBperp,ShockNormalAngle] = thor_QB(tsBSNX,rTHOR,tsB,'perp');
Qpar = QR*QBpar*QV; 
Qperp = QR*QBperp*QV; 
Qpar.data(noCrossing,:)  = 0;
Qperp.data(noCrossing,:) = 0;

%% Turn a 1-min bowshock crossing into a 30-min intervals with corresponding FOM
dtBS = 30*60;
dt = 60; 
tBSparToDownload = 50*30;
tBSperToDownload = 20*30;
% [Qpar30,indQpar30] = thor_30minbs(Qpar,T); % one point is still one minute
% [Qper30,indQper30] = thor_30minbs(Qperp,T);
Qpar30 = thor_q_sort(Qpar,dtBS); % one point is still one minute
Qper30 = thor_q_sort(Qperp,dtBS);

% assign FOMs to sorted lists
FOMpar=Qpar30;FOMpar.data=zeros(size(Qpar30.data));
FOMper=Qper30;FOMper.data=zeros(size(Qper30.data));

tmpQ = Qpar30.data;
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend');
Qlim = [Qlist(tBSparToDownload+1) + (0:-1:-5)*0.1 0.03 0.02 0.01 0];
for iQ = numel(Qlim):-1:1
	tmpFOM = 10-iQ;
	tmpQlim = Qlim(iQ);
	FOMpar.data(tmpQ>tmpQlim)=tmpFOM;
end
tmpQ = Qper30.data;
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend');
Qlim = [Qlist(tBSperToDownload+1) + (0:-1:-5)*0.015 0.03 0.02 0.01 0];
for iQ = numel(Qlim):-1:1
	tmpFOM = 10-iQ;
	tmpQlim = Qlim(iQ);
	FOMper.data(tmpQ>tmpQlim)=tmpFOM;
end
iFOMPerPutToZero=find((FOMpar.data>0) & (FOMper.data > 0));
FOMper.data(iFOMPerPutToZero) = 0;
Qper30.data(iFOMPerPutToZero) = 0;
FOM.data = FOMpar.data + FOMper.data;
iKSR.data(FOMpar.data>0) = 6;
iKSR.data(FOMper.data>0) = 6;

if 0 % plot regions
  %%
isub = 1;
hca = subplot(2,3,isub); isub = isub + 1;
ind=1;
plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'r.')
hca = subplot(2,3,isub); isub = isub + 1;
ind=2;
plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'r.')
hca = subplot(2,3,isub); isub = isub + 1;
ind=3;
plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'r.')
hca = subplot(2,3,isub); isub = isub + 1;
ind=4;
plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'r.')
hca = subplot(2,3,isub); isub = isub + 1;
ind=1; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.'); hold(hca,'on')
ind=2; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.'); 
ind=4; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.');
ind=3; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.'); 
ind=6; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.'); 
hold(hca,'off')

hca = subplot(2,3,isub); isub = isub + 1;
plot(rTHOR.x.data(isWithin15min,:),rTHOR.y.data(isWithin15min,:),'r.'); hold(hca,'on')
ind=6; plot(rTHOR.x.data(find((iKSR.data)==ind),:),rTHOR.y.data(find((iKSR.data)==ind),:),'.'); hold(hca,'off')

end

%% Divide into orbits, and Q-bins
% Quality factor;
% edgesQ = 0.01:0.1:1; % quality factor bins, all Q
% [distQpar30_00,~] = thor_bin(Qpar30,edgesQ,tPerigee);
% [distQper30_00,~] = thor_bin(Qper30,edgesQ,tPerigee);
% edgesQ = 0.8:0.02:1; % quality factor bins, high Q
% [distQpar30_08,~] = thor_bin(Qpar30,edgesQ,tPerigee);
% [distQper30_08,~] = thor_bin(Qper30,edgesQ,tPerigee);
%limQpar1 = sortedQpar(end-40);
%limQpar2 = sortedQpar(end-80);
%limQpar3 = sortedQpar(end-120);

edgesQvar = [1 sortedQpar((1:5)*40)'];
edgesFOM = 9.5:-1:-0.5;
%edgesQvar = [limQpar1 limQpar2 limQpar3]
%edgesQvar = [0.01 0.2 0.4 0.6 0.7 0.8:0.02:1]; % quality factor bins, high Q
% distQpar30_var = thor_bin(Qpar30,edgesQvar,tPerigee);
% distQper30_var = thor_bin(Qper30,edgesQvar,tPerigee);
distFOM = thor_bin(FOM,edgesFOM,tPerigee);

iKSRorbit = thor_bin(iKSR,0.5:1:6.5,tPerigee);

tApogee = iKSRorbit.time;

%% Make plot
if 1
h = irf_plot(3);
  hca = irf_panel('KSR dist total counts');
  tsCumSumQ = irf.ts_scalar(tsDistQpar_resamp30min.time,cumsum(tsDistQpar_resamp30min.data(:,end:-1:1),2));
  hp = irf_patch(hca,tsCumSumQ);
  %irf_plot(hca,tsCumSumKSR);
  hca.YLabel.String = {'Q','#/orbit'};
  labels = arrayfun(@(x,y) {[num2str(x) ' > Q_{||} > ' num2str(y)]}, edgesQ(end:-1:2),edgesQ(end-1:-1:1));
  hpatches=findall(hca,'Type','patch');
  legend(hp,labels,'location','eastoutside')

  hca = irf_panel('qpar dist total counts');
  tsCumSumKSR = irf.ts_scalar(iKSRorbit.time,cumsum(iKSRorbit.data(:,[6 4:-1:1]),2));
  hp = irf_patch(hca,tsCumSumKSR);
  %irf_plot(hca,tsCumSumKSR);
  hca.YLabel.String = {'KSR','#/orbit'};
  legend(hp,{'6-bowshock','4-pristine solar wind','3-foreshock (parker)','2-magnetosheath','1-magnetosphere'},'location','eastoutside')

  hca = irf_panel('KSR + qpar dist total counts');
  tsCumSumKSR_Q = irf.ts_scalar(tsDistKSR_Q.time,cumsum(tsDistKSR_Q.data(:,[14:-1:6 4:-1:1]),2));
  hp = irf_patch(hca,tsCumSumKSR_Q);
  %irf_plot(hca,tsCumSumKSR);
  hca.YLabel.String = {'KSR','#/orbit'};
  %legend(hp,{'6-bowshock','4-pristine solar wind','3-foreshock (parker)','2-magnetosheath','1-magnetosphere'},'location','eastoutside')

  irf_plot_axis_align
end
%% Datavolume
% 1 - magnetosphere, inside magnetopause
% 2 - magnetosheath, outside magnetopause but inside bowshock
% 3 - quasi-parallel foreshock, outside bowshock but inside the
%     region limited by the tangent of IMF B to the bowshock
% 4 - pristine solar wind - outside bowshock and foreshock
% 5 - magnetopause crossing
% 6-end - bowshock crossing: Q: from low to high (6 lowest, end highest)

dt = iKSR.time(2)-iKSR.time(1);

% TM rates with margin
mshTM = 17897*1e-6*1.2; % Gbps
bsTM  = 23572*1e-6*1.2; % Gbps
fsTM  = 16005*1e-6*1.2; % Gbps
pswTM =  6073*1e-6*1.2; % Gbps

% Collected data, sorted by KSR and BS quality factor
dataMSH = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,1)*dt*mshTM);
% distQ   = distQper30_var; dataBSper = irf.ts_scalar(distQ.time,distQ.data(:,:)*dtBS*bsTM);
% distQ   = distQpar30_var; dataBSpar = irf.ts_scalar(distQ.time,distQ.data(:,:)*dtBS*bsTM);
% dataBS  = dataBSpar;
dataBS  = irf.ts_scalar(distFOM.time,distFOM.data(:,:)*dt*bsTM);
dataFS  = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,3)*dt*fsTM);
dataSW  = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,4)*dt*pswTM);

c_eval('ind_phase?  = iKSRorbit.time.tlim(tintPhase?);',1:3);
c_eval('n_phase?    = numel(ind_phase?);'              ,1:3);
c_eval('time_phase? = iKSRorbit.time(ind_phase?);'     ,1:3)

%% What is to be downlinked
% bs (high to low priority) / msh / fs / psw
burstTMperOrbitGbit = 150*0.9; 
phase1_downlink_per_orbit = burstTMperOrbitGbit*[0.8 0.8 0.0 0.0]/sum([0.8 0.8 0.0 0.0]); % Gbit
phase2_downlink_per_orbit = burstTMperOrbitGbit*[0.2 0.2 0.85 0.3]/sum([0.2 0.2 0.85 0.3]); % Gbit
phase3_downlink_per_orbit = burstTMperOrbitGbit*[0.0 0.0 0.15 0.7]/sum([0.0 0.0 0.15 0.7]); % Gbit

tsDownlinkPhase1 = irf.ts_scalar(time_phase1,repmat(phase1_downlink_per_orbit,n_phase1,1)); tsDownlinkPhase1.units = 'Gbit';
tsDownlinkPhase2 = irf.ts_scalar(time_phase2,repmat(phase2_downlink_per_orbit,n_phase2,1)); tsDownlinkPhase2.units = 'Gbit';
tsDownlinkPhase3 = irf.ts_scalar(time_phase3,repmat(phase3_downlink_per_orbit,n_phase3,1)); tsDownlinkPhase3.units = 'Gbit'; 
tsDownlink       = combine(tsDownlinkPhase1,combine(tsDownlinkPhase2,tsDownlinkPhase3));

downlinkMSH = irf.ts_scalar(tsDownlink.time,tsDownlink.data(:,1));
downlinkBS  = irf.ts_scalar(tsDownlink.time,tsDownlink.data(:,2));
downlinkFS  = irf.ts_scalar(tsDownlink.time,tsDownlink.data(:,3));
downlinkSW  = irf.ts_scalar(tsDownlink.time,tsDownlink.data(:,4));


%% OLD Check how much data accumulates onboard the satellite
onboard_memory = 1000; % Gbit
downlink_delay = 1; % data might stay this long before getting downlinked
netDataMSH   = thor_netdata(dataMSH,  downlinkMSH,onboard_memory,downlink_delay);
netDataBS    = thor_netdata(dataBS,   downlinkBS, onboard_memory,downlink_delay);
netDataBSpar = thor_netdata(dataBSpar,downlinkBS, onboard_memory,downlink_delay);
netDataBSper = thor_netdata(dataBSper,downlinkBS, onboard_memory,downlink_delay);
netDataFS    = thor_netdata(dataFS,   downlinkFS, onboard_memory,downlink_delay);
netDataSW    = thor_netdata(dataSW,   downlinkSW, onboard_memory,downlink_delay);
%% Check how much data accumulates onboard the satellite
memorySaved = 1000; % Gbit
tmMSH   = thor_tm(dataMSH,  downlinkMSH,memorySaved);
tmBS    = thor_tm(dataBS,   downlinkBS ,memorySaved);
tmFS    = thor_tm(dataFS,   downlinkFS ,memorySaved);
tmSW    = thor_tm(dataSW,   downlinkSW ,memorySaved);
dataTotal = dataBS;
dataTotal.data = dataBS.data;
dataTotal.data(:,3:10) = dataTotal.data(:,3:10) + dataMSH.data*[0.01 0.01 0.05 0.1 0.1 0.1 0.1 0.53];
dataTotal.data(:,3:10) = dataTotal.data(:,3:10) + dataFS.data*[0.01 0.01 0.05 0.1 0.1 0.1 0.1 0.53];
dataTotal.data(:,3:10) = dataTotal.data(:,3:10) + dataSW.data*[0.025 0.025 0.1 0.1 0.1 0.1 0.1 0.45];
tmTotal = thor_tm(dataTotal,   burstTMperOrbitGbit ,memorySaved);
%% Bowshock datavolume, parallel and perpendicular crossings
h = irf_plot(4);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
 ud.t_start_epoch = dataBS.time(1).epochUnix;
 set(gcf,'userdata',ud);
end

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
	hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,dataTotal.data,1.6,'stacked');
	irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 1000]);
	%labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 9:-1:1);
  legend(hp,labels{:},'location','eastoutside')
	hold(hca,'on')
%   h_downlink = irf_plot(hca,downlinkBS);
%   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
%   hold(hca,'on') 
  hca.YLabel.String = {'New data','Gbit/orbit'};
end
if 1 % downloaded data per orbit
  hca = irf_panel(h,'Onboard data par');
	hp=bar(hca,tmBSpar.downloaded.time.epochUnix-tStartEpoch,tmTotal.downloaded.data,1.5,'stacked');
	grid(hca,'on');
	irf_timeaxis(hca);
	irf_zoom(hca,'y',[0 burstTMperOrbitGbit]);
  labels = arrayfun(@(x) {num2str(x)}, 9:-1:1);
  legend(hp,labels{:},'location','eastoutside')
  hca.YLabel.String = {'Downloaded','Gbit/orbit'};  
end
if 1 % saved data onboard
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,tmBSpar.saved.time.epochUnix-tStartEpoch,tmTotal.saved.data,1.5,'stacked');
  labels = arrayfun(@(x) {num2str(x)}, 9:-1:1);
  irf_zoom(hca,'y',[0 1200]);
	legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
	irf_zoom(hca,'y',[0 memorySaved]);
  hca.YLabel.String = {'Saved onboard','Gbit/orbit'};  
end
if 1 % discarded data onboard
  hca = irf_panel(h,'Discarded data');
	hp=bar(hca,tmBSpar.discarded.time.epochUnix-tStartEpoch,tmTotal.discarded.data,1.5,'stacked');
	irf_zoom(hca,'y',[0 1500]);
	labels = arrayfun(@(x) {num2str(x)}, 9:-1:1);
  legend(hp,labels{:},'location','eastoutside')
  hca.YLabel.String = {'Discarded','Gbit/orbit'};  
end
irf_plot_axis_align

irf_zoom(h,'x',EpochTT([tintPhase1(1).utc; tintPhase3(2).utc]))

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';
%irf_zoom(h,'x',tintTHOR)

%% Write to file
fileID = fopen('/Users/Cecilia/Matlab/irfu-matlab/mission/thor/orbit_coverage/KSR_time_spent.txt','w');
fprintf(fileID,'%6s\n','Time spent in key science regions (minutes)');
fprintf(fileID,'Qpar and Qperp are sorted in bins: %s\n',num2str(edgesQvar));
stringFormat = '%6s %6s %6s %6s %6s %40s %40s';
fprintf(fileID,'%6s %6s %6s %6s %6s %15s %6s %6s %6s %6s %6s %6s %6s %40s\n','orbit','msp','msh','fs','psw','bs: Qpar','bs: Qperp');

numFormat = '%6.0f %6.0f %6.0f %6.0f %6.0f';
for iQ = 1:numel(edgesQvar)
  numFormat = [numFormat ' %6.0f'];
end
for iQ = 1:numel(edgesQvar)
  numFormat = [numFormat ' %6.0f'];
end

printData = [tocolumn(1:iKSRorbit.length) iKSRorbit.data(:,1:4) distQpar30_var.data distQper30_var.data]';
fprintf(fileID,[numFormat '\n'],printData);

%fprintf(fileID,'%6.0f %6.0f\n',[1:tsDistKSR_Q.length;1:tsDistKSR_Q.length]);
fclose(fileID); 

%% Write to excel
csvwrite('/Users/Cecilia/Matlab/irfu-matlab/mission/thor/orbit_coverage/KSR_time_spent',printData')

