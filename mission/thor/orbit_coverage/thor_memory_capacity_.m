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
%datastore('spice','dir','/Users/andris/calc/spice');
units = irf_units;

resampleKernels = 1;
rTHOR_orig = thor_orbit('new2a.bsp',1*3600);
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
% This is how many crossings we want to get down, it should be adjust so
% that we get the appropriate number ofr crossing during the right phase.
% E.g. telemetry is only allocated for BS duing NSP1 and NSP2. Therefore,
% the Qlist below needs to be limieted to these phases, otherwise the
% number of crossing will be less than desired...

% Start by obtaining the times/indices of different phases.
iEndNSP1 = find(rTHOR.abs.data*1e3/units.RE>17,1,'first');
tEndNSP1 = rTHOR.time(iEndNSP1);
iEndNSP2 = find(rTHOR.abs.data*1e3/units.RE>27,1,'first');
tEndNSP2 = rTHOR.time(iEndNSP2);
tNSP12 = rTHOR(1:iEndNSP2).time;
tNSP23 = rTHOR(iEndNSP1:rTHOR.length).time;


% This is the number of required crossings. The categories boundaries will
% be set such that they are obtained during NSP1 & 2
tBSparToDownloadCat1 = 50*30;
tBSperToDownloadCat1 = 20*30;
% Divide the remaining crossings into 3 parts, 20%, 50%, 100%, see below

% [Qpar30,indQpar30] = thor_30minbs(Qpar,T); % one point is still one minute
% [Qper30,indQper30] = thor_30minbs(Qperp,T);
Qpar30 = thor_q_sort(Qpar,dtBS); % one point is still one minute
Qper30 = thor_q_sort(Qperp,dtBS);

% assign FOMs to sorted lists
FOMpar=Qpar30;FOMpar.data=zeros(size(Qpar30.data));
FOMper=Qper30;FOMper.data=zeros(size(Qper30.data));

tmpQ = Qpar30.data; % includes all phases
tmpQ = Qpar30.tlim(tNSP12).data; % limit to NSP1 and 2
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend'); % highest first
Qlist234 = Qlist(tBSparToDownloadCat1+1:end);
%QlimPar = [Qlist(tBSparToDownloadCat1+1) + (0:-1:-5)*0.1 0.03 0.02 0.01 0]; % bins for FOM, assign a given range of Qpar a certain FOM
QlimPar = [Qlist(tBSparToDownloadCat1+1) Qlist234(end*[0.2 0.5])' 0]; % bins for FOM, assign a given range of Qpar a certain FOM
%QlimPar = [Qlist(tBSparToDownloadCat1+1) Qlist(tBSparToDownloadCat2+1) Qlist(tBSparToDownloadCat3+1)]; % bins for FOM, assign a given range of Qpar a certain FOM
QedgesPar = [Qlist(1) QlimPar];
for iQ = numel(QlimPar):-1:1
  tmpFOM = iQ;
  tmpQlim = QlimPar(iQ);
  disp(sprintf('tmpFom=%g, tmpLowerQlimPar = %g',tmpFOM,tmpQlim))
  FOMpar.data(tmpQ>tmpQlim)=tmpFOM;
end

tmpQ = Qper30.data; % includes all phases
%tmpQ = Qper30.tlim(tNSP12).data; % limit to NSP1 and 2
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend');
Qlist234 = Qlist(tBSparToDownloadCat1+1:end);
%QlimPer = [Qlist(tBSperToDownload+1) + (0:-1:-5)*0.015 0.03 0.02 0.01 0];
QlimPer = [Qlist(tBSperToDownloadCat1+1) Qlist234(end*[0.2 0.5])' 0]; % bins for FOM, assign a given range of Qpar a certain FOM
%QlimPer = [Qlist(tBSperToDownloadCat1+1) Qlist(tBSperToDownloadCat2+1) Qlist(tBSperToDownloadCat3+1)]; % bins for FOM, assign a given range of Qpar a certain FOM 0.03 0.02 0.01 0]; % bins for FOM, assign a given range of Qpar a certain FOM
QedgesPer = [Qlist(1) QlimPer];
for iQ = numel(QlimPer):-1:1
  tmpFOM = iQ;
  tmpQlim = QlimPer(iQ);
  disp(sprintf('tmpFom=%g, tmpLowerQlimPer = %g',tmpFOM,tmpQlim))
  FOMper.data(tmpQ>tmpQlim)=tmpFOM;
end
iFOMPerPutToZero=find((FOMpar.data>0) & (FOMper.data > 0));
FOMper.data(iFOMPerPutToZero) = 0;
Qper30.data(iFOMPerPutToZero) = 0;
FOM.data = FOMpar.data + FOMper.data;
iKSR.data(FOMpar.data>0) = 6;
iKSR.data(FOMper.data>0) = 6;

%% Assign FOMs to magnetosheaht time intervals, Q = cosd(theta)^2; theta = angle between THOR and bowshock nose
[Qmsh,angleMSH] = thor_Q_msh(rTHOR); % irf_plot({rTHOR,Qmsh,dataMSH});
FOMmsh = Qmsh;
tmpQ = Qmsh.data;
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend');
QlimMSH = [Qlist(end*[0.25 0.5 0.75])' 0]; % bins for FOM, assign a given range of Qpar a certain FOM
QlimMSH = [0.9 0.8 0.7 0];

for iQ = numel(QlimMSH):-1:1
  tmpFOM = iQ;
  tmpQlim = QlimMSH(iQ);
  disp(sprintf('tmpFom=%g, tmpLowerQlimMSH = %g',tmpFOM,tmpQlim))
  FOMmsh.data(tmpQ>tmpQlim)=tmpFOM;
end
FOMmsh.name = 'FOM MSH';
FOMmsh.data(iKSR.data~=2) = 0; % set non-MSH intervals to FOM 0;

%% Assign FOMs to foreshock...
[Qfs,~] = thor_Q_fs(rTHOR); % irf_plot({rTHOR,Qmsh,dataMSH});
FOMfs = Qfs;
% scatter(rTHOR.x.data(1:20:end,:)*1e3/units.RE,rTHOR.y.data(1:20:end,:)*1e3/units.RE,1,FOMfs.data(1:20:end,:))
tmpQ = Qfs.data;
Qlist = tmpQ(tmpQ>0);
Qlist = sort(Qlist,'descend');
QlimFS = [Qlist(end*[0.25 0.5 0.75])' 0]; % bins for FOM, assign a given range of Qpar a certain FOM
QlimFS = [0.75 0.50 0.25 0];
% scatter(rTHOR.x.data(1:20:end,:)*1e3/units.RE,rTHOR.y.data(1:20:end,:)*1e3/units.RE,1,FOMfs.data(1:20:end,:))
for iQ = numel(QlimFS):-1:1
  tmpFOM = iQ;
  tmpQlim = QlimFS(iQ);
  disp(sprintf('tmpFom=%g, tmpLowerQlimFS = %g',tmpFOM,tmpQlim))
  FOMfs.data(tmpQ>tmpQlim)=tmpFOM;
end
FOMfs.name = 'FOM MSH';
FOMfs.data(iKSR.data~=3) = 0; % set non-FS intervals to FOM 0;

%% plot regions
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

%edgesQvar = [1 sortedQpar((1:5)*40)'];
%edgesFOM = 9.5:-1:-0.5;
%edgesQvar = [limQpar1 limQpar2 limQpar3]
%edgesQvar = [0.01 0.2 0.4 0.6 0.7 0.8:0.02:1]; % quality factor bins, high Q
% distQpar30_var = thor_bin(Qpar30,edgesQvar,tPerigee);
% distQper30_var = thor_bin(Qper30,edgesQvar,tPerigee);
%distFOM = thor_bin(FOM,-0.5:9.5,tPerigee); distFOM.userData = 'FOMSs = [0:9]'; % FOM 9 is the highest quality factor
distFOM = thor_bin(FOM,0.5:1:4.5,tPerigee); distFOM.userData = 'FOMSs = [1:4]'; % FOM 1 is the highest quality factor
distFOM.name = 'FOM BS';

iKSRorbit = thor_bin(iKSR,0.5:1:6.5,tPerigee);
iKSRorbit.userData = iKSR.userData;

distFOM_MSH = thor_bin(FOMmsh,0.5:1:4.5,tPerigee);
distFOM_MSH.name = 'FOM MSH';

distFOM_FS = thor_bin(FOMfs,0.5:1:4.5,tPerigee);
distFOM_FS.name = 'FOM FS';


tApogee = iKSRorbit.time;

iOrbEndNSP2 = find(tApogee<tEndNSP2,1,'last');
iOrbNSP12 = find(tApogee<tEndNSP2);
iOrbEndNSP1 = find(tApogee>tEndNSP1,1,'first');
iOrbNSP23 = find(tApogee>tEndNSP1);
iOrbNSP1 = 1:iOrbEndNSP1;
iOrbNSP2 = (iOrbEndNSP1+1):iOrbEndNSP2;
iOrbNSP3 = (iOrbEndNSP2+1):tApogee.length;

%% Make plot
if 0
  %%
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
% mshTM = 17897*1e-6*1.2; % Gbps
% bsTM  = 23572*1e-6*1.2; % Gbps
% fsTM  = 16005*1e-6*1.2; % Gbps
% pswTM =  6073*1e-6*1.2; % Gbps, for one orbit: pswTMorbit = 6073*1e-6*1.2*(6*24*60*60)

% TM rates with margin, from Payload spreadsheet
mshTM = 14914*1e-6*1.2; % Gbps
bsTM  = 19644*1e-6*1.2; % Gbps
fsTM  = 13337*1e-6*1.2; % Gbps
pswTM =  5061*1e-6*1.2; % Gbps, for one orbit: pswTMorbit = 6073*1e-6*1.2*(6*24*60*60)

% Mission sucess criteria from Payload spreadsheet
mshMSC = 6443;
bsMSC = 2970;
fsMSC = 2881;
swMSC = 3279;

% Collected data, sorted by KSR and BS quality factor
dataMSH = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,2)*dt*mshTM); dataMSH.units = 'Gbit'; % i_2: 'magnetosheath, outside magnetopause but inside bowshock'
dataMSH = irf.ts_scalar(distFOM_MSH.time,distFOM_MSH.data(:,:)*dt*mshTM); dataMSH.units = 'Gbit';
% distQ   = distQper30_var; dataBSper = irf.ts_scalar(distQ.time,distQ.data(:,:)*dtBS*bsTM);
% distQ   = distQpar30_var; dataBSpar = irf.ts_scalar(distQ.time,distQ.data(:,:)*dtBS*bsTM);
% dataBS  = dataBSpar;
dataBS  = irf.ts_scalar(distFOM.time,distFOM.data(:,:)*dt*bsTM); dataBS.units = 'Gbit';% both parallel and perpendicular FOMS combined
%dataFS  = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,3)*dt*fsTM); dataFS.units = 'Gbit';
dataFS  = irf.ts_scalar(distFOM_FS.time,distFOM_FS.data(:,:)*dt*fsTM); dataFS.units = 'Gbit';
dataSW  = irf.ts_scalar(iKSRorbit.time,iKSRorbit.data(:,4)*dt*pswTM); dataSW.units = 'Gbit';

c_eval('ind_phase?  = iKSRorbit.time.tlim(tintPhase?);',1:3);
c_eval('n_phase?    = numel(ind_phase?);'              ,1:3);
c_eval('time_phase? = iKSRorbit.time(ind_phase?);'     ,1:3)

%% What is to be downlinked
% bs (high to low priority) / msh / fs / psw
tmAllocationNSP1 = [0.8 0.8 0.0  0.0];
tmAllocationNSP2 = [0.2 0.2 0.83 0.3];
tmAllocationNSP3 = [0.0 0.0 0.15 0.7];

burstTMperOrbitGbit = 150*0.9;
phase1_downlink_per_orbit = burstTMperOrbitGbit*tmAllocationNSP1/sum(tmAllocationNSP1); % Gbit
phase2_downlink_per_orbit = burstTMperOrbitGbit*tmAllocationNSP2/sum(tmAllocationNSP2); % Gbit
phase3_downlink_per_orbit = burstTMperOrbitGbit*tmAllocationNSP3/sum(tmAllocationNSP3); % Gbit

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

%% OLD Check how much data accumulates onboard the satellite
memorySaved = 12000; % Gbit

tmMSH   = thor_tm(dataMSH, downlinkMSH, memorySaved); tmMSH.collected = dataMSH;
tmBS    = thor_tm(dataBS,  downlinkBS , memorySaved); tmBS.collected = dataBS;
tmFS    = thor_tm(dataFS,  downlinkFS , memorySaved); tFS.collected = dataFS;
tmSW    = thor_tm(dataSW,  downlinkSW , memorySaved); tsSW.collected = dataSW;

% Assign FOM to MSH, FS, and SW data. The FOMs are not assigned to any
% particular time interval, but to a fraction of the data.
% MSH - 120 h in total
% FS - 50 h in total
% SW - 150 h in total
dataTotal = dataBS; % dataBS is the collected datavolume
dataTotal.data = dataBS.data;
dataTotal.data(:,1:4) = dataTotal.data(:,1:4) + dataMSH.data*fomMSH/sum(fomMSH); % grading the importance of other regions data
dataTotal.data(:,1:4) = dataTotal.data(:,1:4) + dataFS.data*fomFS/sum(fomFS);
dataTotal.data(:,1:4) = dataTotal.data(:,1:4) + dataSW.data*fomSW/sum(fomSW);

tmTotal = thor_tm(dataTotal,   burstTMperOrbitGbit ,memorySaved);  tmTotal.collected = irf.ts_scalar(tmTotal.saved.time,dataTotal.data);
%tmTotal = thor_tm_sitl(dataTotal,   burstTMperOrbitGbit ,memorySaved,2);  tmTotal.collected = irf.ts_scalar(tmTotal.saved.time,dataTotal.data);

%% Check how much data accumulates onboard the satellite
memorySaved = 12000; % Gbit

% How much time is in cat 1 msh during nsp1+nsp2
totTmshFOM = sum(distFOM_MSH(iOrbNSP1).data,1); % total time spent in MSH (sorted by FOMs) during NSP1 and NSP2
relTmshFOM = totTmshFOM/sum(totTmshFOM);

totTswFOM = sum(iKSRorbit(iOrbNSP23).data(:,3),1);
relTswFOM = totTswFOM/sum(totTswFOM);

% Assign FOM to MSH, FS, and SW data. The FOMs are not assigned to any
% particular time interval, but to a fraction of the data.
% MSH - 120 h in total
% FS - 50 h in total
% SW - 150 h in total
fomBS = [1.2 1.1 1 1]; fomBS = fomBS/sum(fomBS);
fomMSH = [0.07 0.2 0.3 0.4]; fomMSH = fomMSH/sum(fomMSH);relTmshFOM; % This adjusts the FOM so that the Cat1 is fully acieved during NSP1 and NSP2 % [0.1 0.2 0.3 0.4]; fomMSH = fomMSH/sum(fomMSH); % percentage
fomFS = [0.007 0.05 0.3 0.5]; fomFS = fomFS/sum(fomFS); % percentage
fomSW = [0.03 0.1 0.3 0.4]; fomSW = fomSW/sum(fomSW); % percentage

dataFomMSH = irf.ts_scalar(dataMSH.time,dataMSH.data.*repmat(4*fomMSH,dataMSH.length,1)); % maybe I should weight this with magnetic clock angle
%dataFomFS = irf.ts_scalar(dataFS.time,dataFS.data*fomFS);
dataFomFS = irf.ts_scalar(dataFS.time,dataFS.data.*repmat(4*fomFS,dataFS.length,1));
dataFomSW = irf.ts_scalar(dataSW.time,dataSW.data*fomSW);
dataFomBS = dataBS;
dataFomBS = irf.ts_scalar(dataBS.time,dataBS.data.*repmat(4*fomBS,dataBS.length,1));

if 0
  tmMSH   = thor_tm(dataFomMSH, downlinkMSH, memorySaved); tmMSH.collected = dataFomMSH;
  tmBS    = thor_tm(dataFomBS,  downlinkBS , memorySaved); tmBS.collected = dataFomBS;
  tmFS    = thor_tm(dataFomFS,  downlinkFS , memorySaved); tmFS.collected = dataFomFS;
  tmSW    = thor_tm(dataFomSW,  downlinkSW , memorySaved); tmSW.collected = dataFomSW;
  dataTotal = dataMSH;
  dataTotal.data = dataFomBS.data + dataFomMSH.data + dataFomFS.data + dataFomSW.data;
  tmTotal = thor_tm(dataTotal,   burstTMperOrbitGbit ,memorySaved);  tmTotal.collected = irf.ts_scalar(tmTotal.saved.time,dataTotal.data);
else
  % downlinkDelay = 3;
  % make downlinkDelay into time array, because it varies between the
  % phases, NSP1: 3, NSP2: 2, NSP3: 2

  % make tm buffert also as timeseries: (orbit time - 8 h) * (average telemetry rate)
  % tmAllocationNSP1, tmAllocationNSP2, tmAllocationNSP3
  tRoiNSP1 = (tApogee(iOrbNSP1(end-10))-tApogee(iOrbNSP1(end-11)))-60*60*8; % minus 8 h
  tRoiNSP2 = (tApogee(iOrbNSP2(end-10))-tApogee(iOrbNSP2(end-11)))-60*60*8;
  tRoiNSP3 = (tApogee(iOrbNSP3(end-10))-tApogee(iOrbNSP3(end-11)))-60*60*8;
  atnsp1 = (tmAllocationNSP1(1)*bsTM + tmAllocationNSP1(2)*mshTM + tmAllocationNSP1(3)*fsTM + tmAllocationNSP1(4)*pswTM)/sum(tmAllocationNSP1);
  atnsp2 = (tmAllocationNSP2(1)*bsTM + tmAllocationNSP2(2)*mshTM + tmAllocationNSP2(3)*fsTM + tmAllocationNSP2(4)*pswTM)/sum(tmAllocationNSP2);
  atnsp3 = (tmAllocationNSP3(1)*bsTM + tmAllocationNSP3(2)*mshTM + tmAllocationNSP3(3)*fsTM + tmAllocationNSP3(4)*pswTM)/sum(tmAllocationNSP3);
  avtmNSP1 = atnsp1*tRoiNSP1;
  avtmNSP2 = atnsp2*tRoiNSP2;
  avtmNSP3 = atnsp3*tRoiNSP3;
  nOrbitToStoreData1 = 3;
  nOrbitToStoreData2 = 2;
  nOrbitToStoreData3 = 1;
  c_eval('disp(sprintf(''NSP?: Av TM rate: %.5f Gbps, tOrbit: %5.1f h, Average TM per orbit = %4.0f Gbps, Min # Orbits to store data onboard: %g, Memory buffert needed: %5.0f Gbit'',atnsp?,tRoiNSP?/60/60,avtmNSP?,nOrbitToStoreData?,avtmNSP?*nOrbitToStoreData?))',1:3)
  % bias more towards fs data instead
  %avtmNSP3 = (0*mshTM + 0*bsTM + 0.5*fsTM + 0.5*pswTM)*tRoiNSP3;
  tsDownlinkDelay = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+3;iOrbNSP2'*0+2;iOrbNSP3'*0+1]);
  % how much memory buffert is needed
  %tsAvTM = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+avtmNSP1*nOrbitToStoreData1;iOrbNSP2'*0+avtmNSP2*nOrbitToStoreData2;iOrbNSP3'*0+avtmNSP3*nOrbitToStoreData3]);
  tsAvTM = irf.ts_scalar(dataBS.time,repmat(memorySaved*(1-1/1.75),dataBS.length,1));
  c_eval('relAll? = tmAllocationNSP?/sum(tmAllocationNSP?);',1:3) % bs / msh / fs / psw
  %   tsAvTMmsh = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+avtmNSP1*3*relAll1(2);iOrbNSP2'*0+avtmNSP2*2*relAll2(2);iOrbNSP3'*0+avtmNSP3*1*relAll3(2)]);
  %   tsAvTMbs = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+avtmNSP1*3*relAll1(1);iOrbNSP2'*0+avtmNSP2*2*relAll2(1);iOrbNSP3'*0+avtmNSP3*1*relAll3(1)]);
  %   tsAvTMfs = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+avtmNSP1*3*relAll1(3);iOrbNSP2'*0+avtmNSP2*2*relAll2(3);iOrbNSP3'*0+avtmNSP3*1*relAll3(3)]);
  %   tsAvTMsw  = irf.ts_scalar(dataBS.time,[iOrbNSP1'*0+avtmNSP1*3*relAll1(4);iOrbNSP2'*0+avtmNSP2*2*relAll2(4);iOrbNSP3'*0+avtmNSP3*1*relAll3(4)]);
  %   tmMSH   = thor_tm_sitl___(dataFomMSH, downlinkMSH, memorySaved, tsDownlinkDelay, tsAvTMmsh); tmMSH.collected = dataFomMSH;
  %   tmBS    = thor_tm_sitl___(dataFomBS,  downlinkBS , memorySaved, tsDownlinkDelay, tsAvTMbs); tmBS.collected = dataFomBS;
  %   tmFS    = thor_tm_sitl___(dataFomFS,  downlinkFS , memorySaved, tsDownlinkDelay, tsAvTMfs); tmFS.collected = dataFomFS;
  %   tmSW    = thor_tm_sitl___(dataFomSW,  downlinkSW , memorySaved, tsDownlinkDelay, tsAvTMsw); tmSW.collected = dataFomSW;
  tmMSH   = thor_tm_sitl(dataFomMSH, downlinkMSH, memorySaved, tsDownlinkDelay, tsAvTM); tmMSH.collected = dataFomMSH;
  tmBS    = thor_tm_sitl(dataFomBS,  downlinkBS , memorySaved, tsDownlinkDelay, tsAvTM); tmBS.collected = dataFomBS;
  tmFS    = thor_tm_sitl(dataFomFS,  downlinkFS , memorySaved, tsDownlinkDelay, tsAvTM); tmFS.collected = dataFomFS;
  tmSW    = thor_tm_sitl(dataFomSW,  downlinkSW , memorySaved, tsDownlinkDelay, tsAvTM); tmSW.collected = dataFomSW;

  dataTotal = dataMSH;
  dataTotal.data = dataFomBS.data + dataFomMSH.data + dataFomFS.data + dataFomSW.data;
  tmTotal = thor_tm_sitl(dataTotal,   burstTMperOrbitGbit ,memorySaved, tsDownlinkDelay, tsAvTM);  tmTotal.collected = irf.ts_scalar(tmTotal.saved.time,dataTotal.data);
end

tmBS.saved.units = ''; tmBS.downloaded.units = ''; tmBS.discarded.units = ''; tmBS.collected.units = '';
tmTotal_2.saved = tmMSH.saved + tmBS.saved + tmFS.saved + tmSW.saved;
tmTotal_2.downloaded = tmMSH.downloaded + tmBS.downloaded + tmFS.downloaded + tmSW.downloaded;
tmTotal_2.discarded = tmMSH.discarded + tmBS.discarded + tmFS.discarded + tmSW.discarded;
tmTotal_2.collected = tmMSH.collected + tmBS.collected + tmFS.collected + tmSW.collected;
%tmTotal = thor_tm_sitl(dataTotal,   burstTMperOrbitGbit ,memorySaved,2);  tmTotal.collected = irf.ts_scalar(tmTotal.saved.time,dataTotal.data);

% Test, to see how much tm is downlinked for each region in each phase.
% TM rates with margin
% mshTM = 17897*1e-6*1.2; % Gbps
% bsTM  = 23572*1e-6*1.2; % Gbps
% fsTM  = 16005*1e-6*1.2; % Gbps
% pswTM =  6073*1e-6*1.2; % Gbps, for one orbit: pswTMorbit = 6073*1e-6*1.2*(6*24*60*60)
totMshTM = mshTM*(60*60*120)/1.2;  tsTotMshTM = irf.ts_scalar(tmMSH.downloaded.time([1 end]),totMshTM*[1 1]');
totBsTM = bsTM*(60*60*0.5*70)/1.2; tsTotBsTM = irf.ts_scalar(tmBS.downloaded.time([1 end]),totBsTM*[1 1]');
totFsTM = fsTM*(60*60*50)/1.2;     tsTotFsTM = irf.ts_scalar(tmFS.downloaded.time([1 end]),totFsTM*[1 1]');
totSwTM = pswTM*(60*60*150)/1.2;   tsTotSwTM = irf.ts_scalar(tmSW.downloaded.time([1 end]),totSwTM*[1 1]');
tsTotalTM = tsTotMshTM + tsTotBsTM + tsTotFsTM + tsTotSwTM;

% Common colormap
cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);

%% Total datavolume, FOM, Yellow book, fig 1
h = irf_plot(2);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected Cat1 data per orbit for the different KSRs
  dataCat1 = dataBS;
  dataCat1.data = [tmMSH.collected.data(:,1) tmBS.collected.data(:,1) tmFS.collected.data(:,1) tmSW.collected.data(:,1)];
  hca = irf_panel(h,'Collected data Cat1');
  hp=bar(hca,dataCat1.time.epochUnix-tStartEpoch,dataCat1.data,1.6,'stacked');
  irf_timeaxis(hca);
  %hca.YLabel.String = {'Collected data','Category 1','Gbit/orbit'};
  hca.YLabel.String = {'Collected data','(Gbit/orbit)'};
  %hl = legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside'); %hl.Box = 'off'; hl.Position = [0.1722    0.1084    0.1052    0.1667];
  hl = legend(hca,{'MSH','BS','FS','SW'},'location','northeast'); %hl.Box = 'off'; hl.Position = [0.1722    0.1084    0.1052    0.1667];
end
if 0 % Saved data per orbit
  hca = irf_panel(h,'Saved data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.saved.data,2.5,'stacked');
  hold(hca,'on')
  irf_plot(hca,tsAvTM*-1+memorySaved,'k');
  %tsPatch = irf.ts_scalar(tsAvTM)
  hpatch = irf_patch(hca,{tsAvTM*-1+memorySaved,memorySaved}); hpatch.FaceColor = [0.8 0.8 0.8]; hpatch.EdgeColor = 0.2*[1 1 1];
  hold(hca,'off')
  hl_ = irf_legend(hca,{'Memory allocated for'},[0.01, 0.85],'k'); hl_.Color = 0.3*[1 1 1];
  hl_ = irf_legend(hca,{'newly collected data'},[0.01, 0.7],'k'); hl_.Color = 0.3*[1 1 1];

  %irf_legend(hca,{'Data with assigned FOM'},[0.01, 0.55],'k')
  %irf_legend(hca,{'placed in download queue'},[0.01, 0.3],'k')

  %ht = irf_legend(hca,{'Data with assigned FOM'},[0.98, 0.3],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  %ht = irf_legend(hca,{'placed in download queue'},[0.98, 0.1],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat. ' num2str(x)]}, 1:4);
  %labels = {labels{:},'m'}
  legend([hp],labels{:},'location','eastoutside')


  hca.YLabel.String = {'Saved','onboard','Gbit/orbit'};
  hca.YLabel.String = {'Saved data in','download queue','Gbit'};
  hca.YLabel.String = {'Onboard data','Gbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded Cat1 data');
  set(hca,'ColorOrder',cmap)
  tmTotalCat1 = irf.ts_scalar(tmTotal.downloaded.time,[tmMSH.downloaded.data(:,1), tmBS.downloaded.data(:,1), tmFS.downloaded.data(:,1) tmSW.downloaded.data(:,1)]);
  hh=irf_plot(hca,tmTotalCat1.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  if 0
    iMisReqMSH = find(cumsum(tmMSH.downloaded.data(:,1),1)>tsTotMshTM.data(1,1),1,'first'); if isempty(iMisReqMSH); iMisReqMSH = tmMSH.downloaded.length; end
    iMisReqBS = find(cumsum(tmBS.downloaded.data(:,1),1)>tsTotBsTM.data(1,1),1,'first'); if isempty(iMisReqBS); iMisReqBS = tmMSH.downloaded.length; end
    iMisReqFS = find(cumsum(tmFS.downloaded.data(:,1),1)>tsTotFsTM.data(1,1),1,'first'); if isempty(iMisReqFS); iMisReqFS = tmMSH.downloaded.length; end
    iMisReqSW = find(cumsum(tmSW.downloaded.data(:,1),1)>tsTotSwTM.data(1,1),1,'first'); if isempty(iMisReqSW); iMisReqSW = tmMSH.downloaded.length; end
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqMSH),tsTotMshTM.data(1,1)),'*','color',cmap(1,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqBS),tsTotBsTM.data(1,1)),'*','color',cmap(2,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqFS),tsTotFsTM.data(1,1)),'*','color',cmap(3,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqSW),tsTotSwTM.data(1,1)),'*','color',cmap(4,:));
  else
    iMisReqMSH = find(cumsum(tmMSH.downloaded.data(:,1),1)>mshMSC,1,'first'); if isempty(iMisReqMSH); iMisReqMSH = tmMSH.downloaded.length; end
    iMisReqBS = find(cumsum(tmBS.downloaded.data(:,1),1)>bsMSC,1,'first'); if isempty(iMisReqBS); iMisReqBS = tmMSH.downloaded.length; end
    iMisReqFS = find(cumsum(tmFS.downloaded.data(:,1),1)>fsMSC,1,'first'); if isempty(iMisReqFS); iMisReqFS = tmMSH.downloaded.length; end
    iMisReqSW = find(cumsum(tmSW.downloaded.data(:,1),1)>swMSC,1,'first'); if isempty(iMisReqSW); iMisReqSW = tmMSH.downloaded.length; end
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqMSH),mshMSC),'*','color',cmap(1,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqBS),bsMSC),'*','color',cmap(2,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqFS),fsMSC),'*','color',cmap(3,:));
    irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqSW),swMSC),'*','color',cmap(4,:));
  end
  %   irf_plot(hca,tsTotMshTM,'--','color',cmap(1,:));
  %   irf_plot(hca,tsTotBsTM,'--','color',cmap(2,:));
  %   irf_plot(hca,tsTotFsTM,'--','color',cmap(3,:));
  %   irf_plot(hca,tsTotSwTM,'--','color',cmap(4,:));
  hold(hca,'off')
  %hca.YLabel.String = {'Accumulated','downloaded data','Category 1','Gbit'};
  hca.YLabel.String = {'Accumulated','downloaded data','(Gbit)'};
  %legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside')
  %hl = legend(hca,{'MSH','BS','FS','SW'},'location','west');  hl.Position = [0.1722    0.1084    0.1052    0.1667]; hl.Box = 'off';
  %irf_legend(hca,{'data required'},[0.05, 0.96],'k')
  irf_legend(hca,{'Mission success criteria'},[0.05, 0.96],'k')
  hstar = irf_legend(hca,{'*'},[0.01, 0.95],'k'); hstar.FontSize = 20;
end

c_eval('irf_zoom(h(?),''y'',[0 8000]);',2)
set(h(1),'YLim',[0 400]);
cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);
c_eval('colormap(h(?),cmap);',[1:2])
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''on'';',[1:2])

linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)
c_eval('h(?).Position(4) = h(?).Position(4)*0.9;',1:2)
h(1).Position(2) = h(1).Position(2)*0.93;
%h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';
h(1).Title.String = 'THOR - FOM and Burst TM datavolumes, Cat. 1';

%% Total datavolume, FOM, Yellow book, fig 2
h = irf_plot(1);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Onboard data
  hca = irf_panel(h,'Saved data');
  hp = bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.saved.data,2.5,'stacked');
  hold(hca,'on')
  irf_plot(hca,tsAvTM*-1+memorySaved,'k');
  % Patch showing what memory data is unavailable
  %hpatch = irf_patch(hca,{tsAvTM*-1+memorySaved,memorySaved}); hpatch.FaceColor = [0.8 0.8 0.8]; hpatch.EdgeColor = 0.2*[1 1 1];
  % Patch showing what available memory data is margin
  hpatch2 = irf_patch(hca,{tsAvTM*-1+memorySaved,memorySaved*(1-1/1.75)*0.25});
  hpatch2.FaceColor = [1 1 1];
  hpatch2.EdgeColor = 0.0*[1 1 1];
  hpatch2.FaceAlpha = 0.3;

  %hl_ = irf_legend(hca,{'Memory allocated for'},[0.01, 0.95],'k'); hl_.Color = 0.3*[1 1 1];
  %hl_ = irf_legend(hca,{'newly collected data'},[0.01, 0.8],'k'); hl_.Color = 0.3*[1 1 1];
  if 1
    hl_ = irf_legend(hca,{'Memory allocated for newly collected data'},[0.95, 0.8],'k'); hl_.Color = 0.3*[1 1 1];
    arrow([60*60*24*4*365 3000],[60*60*24*4*365 memorySaved*(1-1/1.75)*0.25])
    arrow([60*60*24*4*365 4500],[60*60*24*4*365 memorySaved*(1/1.75)])
    hl_ = irf_legend(hca,{'Margin'},[0.95, 0.3],'k'); hl_.Color = 0.3*[1 1 1];
    arrow([60*60*24*4*365 8500],[60*60*24*4*365 memorySaved*(1/1.75)])
    arrow([60*60*24*4*365 10000],[60*60*24*4*365 memorySaved])
    labels = arrayfun(@(x) {['Cat. ' num2str(x)]}, 1:4);
    %labels = arrayfun(@(x) {['FOM ' num2str(x)]}, 1:4);
    legend([hp],labels{:},'location','northwest')
  else
    hl_ = irf_legend(hca,{'Memory allocated for newly collected data'},[0.01, 0.8],'k'); hl_.Color = 0.3*[1 1 1];
    arrow([60*60*24*30*2 3000],[60*60*24*30*2 memorySaved*(1-1/1.75)*0.25])
    arrow([60*60*24*30*2 4500],[60*60*24*30*2 memorySaved*(1/1.75)])
    hl_ = irf_legend(hca,{'Margin'},[0.01, 0.3],'k'); hl_.Color = 0.3*[1 1 1];
    arrow([60*60*24*30*2 8500],[60*60*24*30*2 memorySaved*(1/1.75)])
    arrow([60*60*24*30*2 10000],[60*60*24*30*2 memorySaved])
    labels = arrayfun(@(x) {['Cat. ' num2str(x)]}, 1:4);
    legend([hp],labels{:},'location','northeast')
  end
  hold(hca,'off')
  %irf_legend(hca,{'Data with assigned FOM'},[0.01, 0.55],'k')
  %irf_legend(hca,{'placed in download queue'},[0.01, 0.3],'k')
  %ht = irf_legend(hca,{'Data with assigned FOM'},[0.98, 0.3],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  %ht = irf_legend(hca,{'placed in download queue'},[0.98, 0.1],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  irf_timeaxis(hca);


  hca.YLabel.String = {'Saved','onboard','Gbit/orbit'};
  hca.YLabel.String = {'Saved data in','download queue','Gbit'};
  hca.YLabel.String = {'Memory','(Gbit)'};
  hca.YLim = [0 12000];
end


cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);
c_eval('colormap(h(?),cmap);',[1])
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',[1])
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''on'';',[1])

h(1).Position(4) = h(1).Position(4)*0.9;
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';

%% Total datavolume, FOM, Yellow book suggestion
h = irf_plot(3);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected Cat1 data per orbit for the different KSRs
  dataCat1 = dataBS;
  dataCat1.data = [tmMSH.collected.data(:,1) tmBS.collected.data(:,1) tmFS.collected.data(:,1) tmSW.collected.data(:,1)];
  hca = irf_panel(h,'Collected data Cat1');
  hp=bar(hca,dataCat1.time.epochUnix-tStartEpoch,dataCat1.data,2.5,'stacked');
  irf_timeaxis(hca);
  hca.YLabel.String = {'Collected data','Category 1','Gbit/orbit'};
  hl = legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside'); %hl.Box = 'off'; hl.Position = [0.1722    0.1084    0.1052    0.1667];
end
if 0 % Collected data per orbit
  hca = irf_panel(h,'Collected data MSh');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmMSH.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Coll. MSH','Gbit/orbit'};
end
if 0 % Collected data per orbit
  hca = irf_panel(h,'Collected data BS');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmBS.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Coll. BS','Gbit/orbit'};
end
if 0 % Collected data per orbit
  hca = irf_panel(h,'Collected data FS');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmFS.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Coll. FS','Gbit/orbit'};
end
if 0 % Collected data per orbit
  hca = irf_panel(h,'Collected data SW');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmSW.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Coll. SW','Gbit/orbit'};
end
if 0 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Collected','Gbit/orbit'};
end
if 0 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Discarded','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.saved.data,2.5,'stacked');
  hold(hca,'on')
  irf_plot(hca,tsAvTM*-1+memorySaved,'k');
  %tsPatch = irf.ts_scalar(tsAvTM)
  hpatch = irf_patch(hca,{tsAvTM*-1+memorySaved,memorySaved}); hpatch.FaceColor = [0.8 0.8 0.8]; hpatch.EdgeColor = 0.2*[1 1 1];
  hold(hca,'off')
  hl_ = irf_legend(hca,{'Memory allocated for'},[0.01, 0.85],'k'); hl_.Color = 0.3*[1 1 1];
  hl_ = irf_legend(hca,{'newly collected data'},[0.01, 0.7],'k'); hl_.Color = 0.3*[1 1 1];

  %irf_legend(hca,{'Data with assigned FOM'},[0.01, 0.55],'k')
  %irf_legend(hca,{'placed in download queue'},[0.01, 0.3],'k')

  %ht = irf_legend(hca,{'Data with assigned FOM'},[0.98, 0.3],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  %ht = irf_legend(hca,{'placed in download queue'},[0.98, 0.1],'w'); ht.Color = [1 1 1]; ht.FontSize = 11;
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat. ' num2str(x)]}, 1:4);
  %labels = {labels{:},'m'}
  legend([hp],labels{:},'location','eastoutside')


  hca.YLabel.String = {'Saved','onboard','Gbit/orbit'};
  hca.YLabel.String = {'Saved data in','download queue','Gbit'};
  hca.YLabel.String = {'Onboard data','Gbit'};
end
if 0 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  hca.YLabel.String = {'Downloaded','Gbit/orbit'};
end
if 0
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmTotal.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotalTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 1
  hca = irf_panel(h,'Accumulated downloaded Cat1 data');
  set(hca,'ColorOrder',cmap)
  tmTotalCat1 = irf.ts_scalar(tmTotal.downloaded.time,[tmMSH.downloaded.data(:,1), tmBS.downloaded.data(:,1), tmFS.downloaded.data(:,1) tmSW.downloaded.data(:,1)]);
  hh=irf_plot(hca,tmTotalCat1.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  iMisReqMSH = find(cumsum(tmMSH.downloaded.data(:,1),1)>tsTotMshTM.data(1,1),1,'first'); if isempty(iMisReqMSH); iMisReqMSH = tmMSH.downloaded.length; end
  iMisReqBS = find(cumsum(tmBS.downloaded.data(:,1),1)>tsTotBsTM.data(1,1),1,'first'); if isempty(iMisReqBS); iMisReqBS = tmMSH.downloaded.length; end
  iMisReqFS = find(cumsum(tmFS.downloaded.data(:,1),1)>tsTotFsTM.data(1,1),1,'first'); if isempty(iMisReqFS); iMisReqFS = tmMSH.downloaded.length; end
  iMisReqSW = find(cumsum(tmSW.downloaded.data(:,1),1)>tsTotSwTM.data(1,1),1,'first'); if isempty(iMisReqSW); iMisReqSW = tmMSH.downloaded.length; end
  %   irf_plot(hca,tsTotMshTM,'--','color',cmap(1,:));
  %   irf_plot(hca,tsTotBsTM,'--','color',cmap(2,:));
  %   irf_plot(hca,tsTotFsTM,'--','color',cmap(3,:));
  %   irf_plot(hca,tsTotSwTM,'--','color',cmap(4,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqMSH),tsTotMshTM.data(1,1)),'*','color',cmap(1,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqBS),tsTotBsTM.data(1,1)),'*','color',cmap(2,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqFS),tsTotFsTM.data(1,1)),'*','color',cmap(3,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqSW),tsTotSwTM.data(1,1)),'*','color',cmap(4,:));
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Category 1','Gbit'};
  legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside')
  %hl = legend(hca,{'MSH','BS','FS','SW'},'location','west');  hl.Position = [0.1722    0.1084    0.1052    0.1667]; hl.Box = 'off';
  irf_legend(hca,{'data required'},[0.05, 0.96],'k')
  hstar = irf_legend(hca,{'*'},[0.01, 0.95],'k'); hstar.FontSize = 20;
end
if 0
  hca = irf_panel(h,'Thor orbit');
  irf_plot(hca,rTHOR.resample(rTHOR.time(1:5:end)))
end
if 0
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end

c_eval('irf_zoom(h(?),''y'',[0 7000]);',2)
set(h(1),'YLim',[0 400]);
cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);
c_eval('colormap(h(?),cmap);',[1:2])
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',[1:3])
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''on'';',[1:3])

linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';

%% Total datavolume, FOM
h = irf_plot(7);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Collected','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Discarded','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.saved.data,1.6,'stacked');
  irf_timeaxis(hca);labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Saved','onboard','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  hca.YLabel.String = {'Downloaded','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmTotal.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotalTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 1
  hca = irf_panel(h,'Accumulated downloaded Cat1 data');
  set(hca,'ColorOrder',cmap)
  tmTotalCat1 = irf.ts_scalar(tmTotal.downloaded.time,[tmMSH.downloaded.data(:,1), tmBS.downloaded.data(:,1), tmFS.downloaded.data(:,1) tmSW.downloaded.data(:,1)]);
  hh=irf_plot(hca,tmTotalCat1.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  iMisReqMSH = find(cumsum(tmMSH.downloaded.data(:,1),1)>tsTotMshTM.data(1,1),1,'first'); if isempty(iMisReqMSH); iMisReqMSH = tmMSH.downloaded.length; end
  iMisReqBS = find(cumsum(tmBS.downloaded.data(:,1),1)>tsTotBsTM.data(1,1),1,'first'); if isempty(iMisReqBS); iMisReqBS = tmMSH.downloaded.length; end
  iMisReqFS = find(cumsum(tmFS.downloaded.data(:,1),1)>tsTotFsTM.data(1,1),1,'first'); if isempty(iMisReqFS); iMisReqFS = tmMSH.downloaded.length; end
  iMisReqSW = find(cumsum(tmSW.downloaded.data(:,1),1)>tsTotSwTM.data(1,1),1,'first'); if isempty(iMisReqSW); iMisReqSW = tmMSH.downloaded.length; end
  %   irf_plot(hca,tsTotMshTM,'--','color',cmap(1,:));
  %   irf_plot(hca,tsTotBsTM,'--','color',cmap(2,:));
  %   irf_plot(hca,tsTotFsTM,'--','color',cmap(3,:));
  %   irf_plot(hca,tsTotSwTM,'--','color',cmap(4,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqMSH),tsTotMshTM.data(1,1)),'*','color',cmap(1,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqBS),tsTotBsTM.data(1,1)),'*','color',cmap(2,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqFS),tsTotFsTM.data(1,1)),'*','color',cmap(3,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqSW),tsTotSwTM.data(1,1)),'*','color',cmap(4,:));
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Cat 1','Gbit'};
  legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside')
end
if 0
  hca = irf_panel(h,'Thor orbit');
  irf_plot(hca,rTHOR.resample(rTHOR.time(1:5:end)))
end
if 1
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end
%

c_eval('irf_zoom(h(?),''y'',[0 6000]);',3)
cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);
c_eval('colormap(h(?),cmap);',[1 2 3 4])

h(1).YScale = 'lin';
linkaxes(h,'x')
irf_plot_axis_align
%irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';

%% Total_2 datavolume, FOM
h = irf_plot(7);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal_2.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Collected','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal_2.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Discarded','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal_2.saved.data,1.6,'stacked');
  irf_timeaxis(hca);
  %irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Saved onboard','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmTotal_2.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  %hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Downloaded','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmTotal_2.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotalTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 1
  hca = irf_panel(h,'Accumulated downloaded Cat1 data');
  set(hca,'ColorOrder',cmap)
  tmTotalCat1 = irf.ts_scalar(tmTotal.downloaded.time,[tmMSH.downloaded.data(:,1), tmBS.downloaded.data(:,1), tmFS.downloaded.data(:,1) tmSW.downloaded.data(:,1)]);
  hh=irf_plot(hca,tmTotalCat1.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  iMisReqMSH = find(cumsum(tmMSH.downloaded.data(:,1),1)>tsTotMshTM.data(1,1),1,'first');
  iMisReqBS = find(cumsum(tmBS.downloaded.data(:,1),1)>tsTotBsTM.data(1,1),1,'first');
  iMisReqFS = find(cumsum(tmFS.downloaded.data(:,1),1)>tsTotFsTM.data(1,1),1,'first');
  iMisReqSW = find(cumsum(tmSW.downloaded.data(:,1),1)>tsTotSwTM.data(1,1),1,'first');
  %   irf_plot(hca,tsTotMshTM,'--','color',cmap(1,:));
  %   irf_plot(hca,tsTotBsTM,'--','color',cmap(2,:));
  %   irf_plot(hca,tsTotFsTM,'--','color',cmap(3,:));
  %   irf_plot(hca,tsTotSwTM,'--','color',cmap(4,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqMSH),tsTotMshTM.data(1,1)),'*','color',cmap(1,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqBS),tsTotBsTM.data(1,1)),'*','color',cmap(2,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqFS),tsTotFsTM.data(1,1)),'*','color',cmap(3,:));
  irf_plot(hca,irf.ts_scalar(tmMSH.downloaded.time(iMisReqSW),tsTotSwTM.data(1,1)),'*','color',cmap(4,:));
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded, Cat 1','Gbit'};
  legend(hca,{'MSH','BS','FS','SW'},'location','eastoutside')
end
if 1
  hca = irf_panel(h,'Thor orbit');
  irf_plot(hca,rTHOR.resample(rTHOR.time(1:5:end)))
end

%
for iP = [3]
  %irf_zoom(h(iP),'y',[0 1001]);
end
cmap = mms_colors('2b34'); cmap(2,:) = [1 0.8 0]; %cmap = cmap(end:-1:1,:);
%cmap = colormap('hot');
for iP = [1 2 3 4]
  colormap(h(iP),cmap);
end

h(1).YScale = 'lin';
linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes';

%% MSH datavolume, FOM
h = irf_plot(6);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmMSH.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {['Cat ' num2str(x)]}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Collected','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmMSH.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 12000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Discarded','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmMSH.saved.data,1.6,'stacked');
  irf_timeaxis(hca);
  %irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Saved onboard','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmMSH.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  %irf_zoom(hca,'y',[0 120]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Downloaded','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmMSH.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotMshTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 1
  %%
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000*dt];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end

c_eval('colormap(h(?),cmap);',1:4)
c_eval('h(?).YLim = [0 1000];',[1 2])
c_eval('h(?).YLim = [0 13000];',[3])

h(1).YScale = 'lin';
linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes for MSH';

%% Bowshock datavolume, parallel and perpendicular crossings
h = irf_plot(6);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmBS.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 1000]);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'New data','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Collected data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmBS.saved.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Saved data','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmBS.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Downloaded data','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmBS.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 500]);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Discarded data','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmBS.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotBsTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 0
  hca = irf_panel(h,'Key Science region index');
  hh = irf_plot(hca,{iKSR(find(iKSR.data==1)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==2)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==3)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==4)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==6)).resample(iKSR.time(1:10:end))},'comp','none');
  c_eval('hh.Children(?).Marker =''.'';',1:5)
  hca.YLabel.String = 'KSR index';
  legend(hca,{'MSP','MSH','FS','PSW','BS'},'location','eastoutside')
  hca.YLim = [0 7];
end
if 1
  %%
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000*dt];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end



c_eval('colormap(h(?),cmap);',1:4)
c_eval('h(?).YLim = [0 500];',[1 4])
h(2).YLim = [0 7000];
linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes, Bowshock';

%% Foreshock datavolume
h = irf_plot(6);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmFS.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'New data','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmFS.saved.data,1.6,'stacked');
  irf_timeaxis(hca);
  %irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Saved data','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmFS.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);
  %irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'Downloaded data','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmFS.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Discarded data','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmFS.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotFsTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 0
  hca = irf_panel(h,'Key Science region index');
  hh = irf_plot(hca,{iKSR(find(iKSR.data==1)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==2)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==3)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==4)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==6)).resample(iKSR.time(1:10:end))},'comp','none');
  c_eval('hh.Children(?).Marker =''.'';',1:5)
  hca.YLabel.String = 'KSR index';
  legend(hca,{'MSP','MSH','FS','PSW','BS'},'location','eastoutside')
  hca.YLim = [0 7];
end
if 1
  %%
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000*dt];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end

c_eval('h(?).YLim = [0 8000];',1:2)
%c_eval('h(?).YLim = [0 13000];',1:2)
for iP = [1 2 3 4]
  colormap(h(iP),cmap);
end
linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes, Foreshock';

%% Solar wind datavolume,
h = irf_plot(6);
ud = get(gcf,'userdata');
if ~isfield(ud,'t_start_epoch')
  ud.t_start_epoch = dataBS.time(1).epochUnix;
  set(gcf,'userdata',ud);
end
tStartEpoch = ud.t_start_epoch;

if 1 % Collected data per orbit
  hca = irf_panel(h,'Collected data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmSW.collected.data,1.6,'stacked');
  irf_timeaxis(hca);
  irf_zoom(hca,'y',[0 1000]);
  %labels = arrayfun(@(x,y) {[num2str(x,2) ' > Q_{||} > ' num2str(y,2)]}, edgesQvar(1:1:end-1),edgesQvar(2:1:end));
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  %   h_downlink = irf_plot(hca,downlinkBS);
  %   irf_legend(hca,{'- downlink/orbit'},[0.02 0.95],'k')
  %   hold(hca,'on')
  hca.YLabel.String = {'New data','Gbit/orbit'};
end
if 1 % Saved data per orbit
  hca = irf_panel(h,'Saved data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmSW.saved.data,1.6,'stacked');
  irf_timeaxis(hca);labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Saved data','Gbit/orbit'};
end
if 1 % Downloaded data per orbit
  hca = irf_panel(h,'Downloaded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmSW.downloaded.data,1.6,'stacked');
  irf_timeaxis(hca);labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Downloaded data','Gbit/orbit'};
end
if 1 % Discarded data per orbit
  hca = irf_panel(h,'Discarded data par');
  hp=bar(hca,dataBS.time.epochUnix-tStartEpoch,tmSW.discarded.data,1.6,'stacked');
  irf_timeaxis(hca);
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  legend(hp,labels{:},'location','eastoutside')
  hold(hca,'on')
  hca.YLabel.String = {'Discarded data','Gbit/orbit'};
end
if 1
  hca = irf_panel(h,'Accumulated downloaded data');
  set(hca,'ColorOrder',cmap)
  hh=irf_plot(hca,tmSW.downloaded.cumsum('t'));
  c_eval('hh(?).Color = cmap(?,:);',1:4);
  hold(hca,'on')
  irf_plot(hca,tsTotSwTM,'--k');
  hold(hca,'off')
  hca.YLabel.String = {'Accumulated','downloaded','Gbit'};
  labels = arrayfun(@(x) {num2str(x)}, 1:4);
  labels{5} = 'mis. req.';
  legend(hca,labels{:},'location','eastoutside')
end
if 0
  hca = irf_panel(h,'Key Science region index');
  hh = irf_plot(hca,{iKSR(find(iKSR.data==1)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==2)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==3)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==4)).resample(iKSR.time(1:10:end)),...
    iKSR(find(iKSR.data==6)).resample(iKSR.time(1:10:end))},'comp','none');
  c_eval('hh.Children(?).Marker =''.'';',1:5)
  hca.YLabel.String = 'KSR index';
  legend(hca,{'MSP','MSH','FS','PSW','BS'},'location','eastoutside')
  hca.YLim = [0 7];
end
if 1
  hca = irf_panel(h,'Region index');
  cmap_ = [0.9 0.9 0.9; mms_colors('2'); mms_colors('3'); mms_colors('4'); [0 0 0]; mms_colors('b')];
  hp = irf_patch(hca,iKSRorbit.cumsum(2)*dt/60/60);
  legend(hp,{'MSP','MSH','FS','PSW','MP','BS'},'location','eastoutside')
  c_eval('hp(?).FaceColor = cmap_(?,:); hp(?).FaceAlpha = 1;',1:6)
  %hca.YLim = [0 13000*dt];
  hca.YLabel.String = {'KSR','dwell time','hours'};
end

c_eval('h(?).YLim = [0 7000];',1:2)
for iP = [1 2 3 4]
  colormap(h(iP),cmap);
end
linkaxes(h,'x')
irf_plot_axis_align
irf_zoom(h,'x',tmTotal.saved.time)

h(1).Title.String = 'THOR - FOM and Burst TM datavolumes, Solar wind';

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

