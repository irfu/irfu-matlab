% Calculate and plot electron and ion pitch angle distributions from L1b particle 
% brst data. 
% Written by D. B. Graham.

ic = 3; % Spacecraft number

tint = irf.tint('2015-12-30T00:29:44.00Z/2015-12-30T00:30:34.00Z');

%% Load Data
if 1,
tic;
c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ic);
c_eval('energye0=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy0'',tint);',ic);
c_eval('energye1=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy1'',tint);',ic);
c_eval('phie=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_phi'',tint);',ic);
c_eval('thetae=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_theta'',tint);',ic);
c_eval('stepTablee=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_stepTable_parity'',tint);',ic);
toc;
tic;
c_eval('disti = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint);',ic);
c_eval('energyi0=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_energy0'',tint);',ic);
c_eval('energyi1=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_energy1'',tint);',ic);
c_eval('phii=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_phi'',tint);',ic);
c_eval('thetai=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_theta'',tint);',ic);
c_eval('stepTablei=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_stepTable_parity'',tint);',ic);
toc;
end

%c_eval('tmpDataObj = dataobj(''/Volumes/DansHD2/data/mms?/fpi/brst/l1b/des-dist/2015/12/30/mms?_fpi_brst_l1b_des-dist_20151230002944_v1.1.0.cdf'');',ic);
%c_eval('diste = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_brstSkyMap_dist''));',ic);
%c_eval('energye0 = get_variable(tmpDataObj,''mms?_des_brstSkyMap_energy0'');',ic);
%c_eval('energye1 = get_variable(tmpDataObj,''mms?_des_brstSkyMap_energy1'');',ic);
%c_eval('stepTablee = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_stepTable_parity''));',ic);
%c_eval('phie = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_brstSkyMap_phi''));',ic);
%c_eval('thetae = get_variable(tmpDataObj,''mms?_des_brstSkyMap_theta'');',ic);

%c_eval('tmpDataObj = dataobj(''/Volumes/DansHD2/data/mms?/fpi/brst/l1b/dis-dist/2015/12/30/mms?_fpi_brst_l1b_dis-dist_20151230002944_v1.1.0.cdf'');',ic);
%c_eval('disti = mms.variable2ts(get_variable(tmpDataObj,''mms?_dis_brstSkyMap_dist''));',ic);
%c_eval('energyi0 = get_variable(tmpDataObj,''mms?_dis_brstSkyMap_energy0'');',ic);
%c_eval('energyi1 = get_variable(tmpDataObj,''mms?_dis_brstSkyMap_energy1'');',ic);
%c_eval('stepTablei = mms.variable2ts(get_variable(tmpDataObj,''mms?_dis_stepTable_parity''));',ic);
%c_eval('phii = mms.variable2ts(get_variable(tmpDataObj,''mms?_dis_brstSkyMap_phi''));',ic);
%c_eval('thetai = get_variable(tmpDataObj,''mms?_dis_brstSkyMap_theta'');',ic);
%tint = irf.tint(diste.time.start.utc,diste.time.stop.utc);

c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_l2pre'',''mms?_dfg_srvy_l2pre_dmpa'',tint);',ic);

c_eval('SCpot=mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_psp'',tint);',ic);
offset1 = 1.3; offset2 = 1.5; offset3 = 1.2; offset4 = 0.0; %For v1 data
c_eval('SCpot.data = -SCpot.data*1.2+offset?;',ic);

energye0 = energye0.data;
energye1 = energye1.data;
energyi0 = energyi0.data;
energyi1 = energyi1.data;
thetae = thetae.data;
thetai = thetai.data;

%% Compute moments
Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
mp = Units.mp;
imoments = mms.psd_moments(disti,phii,thetai,stepTablei,energyi0,energyi1,SCpot,'ion');
ni = imoments.n_psd;
Vi = imoments.V_psd;
Ti = imoments.T_psd;
Tipp = mms.rotate_tensor(Ti,'fac',Bxyz,'pp'); 
Tiparperp = TSeries(Ti.time,[Tipp.data(:,1,1) Tipp.data(:,2,2)]);

emoments = mms.psd_moments(diste,phie,thetae,stepTablee,energye0,energye1,SCpot,'electron');
ne = emoments.n_psd;
Ve = emoments.V_psd;
Te = emoments.T_psd;
Tepp = mms.rotate_tensor(Te,'fac',Bxyz,'pp'); 
Teparperp = TSeries(Te.time,[Tepp.data(:,1,1) Tepp.data(:,2,2)]);

diste.data = diste.data*1e30; % Unit conversion
disti.data = disti.data*1e30;

energyspec = ones(length(diste.time),1)*energye0;
for ii = 1:length(diste.time);
    if stepTablee.data(ii),
        energyspec(ii,:) = energye1;
    end
end

energyspeci = ones(length(disti.time),1)*energyi0;
for ii = 1:length(disti.time);
    if stepTablei.data(ii),
        energyspeci(ii,:) = energyi1;
    end
end

% define angles
dangle = pi/16;
lengthphi = 32;

z2 = ones(lengthphi,1)*sind(thetae);
solida = dangle*dangle*z2;
allsolidi = zeros(size(disti.data));
allsolide = zeros(size(diste.data));

for ii = 1:length(disti.time);
    for jj=1:length(energyi0);
        allsolidi(ii,jj,:,:) = solida;
    end
end

for ii = 1:length(diste.time);
    for jj=1:length(energye0);
        allsolide(ii,jj,:,:) = solida;
    end
end

distis = disti.data.*allsolidi;
distes = diste.data.*allsolide;

% Electron analysis - OMNI
PSDomni = zeros(length(diste.time),length(energye0));
for ii = 1:length(diste.time);
    disttemp = squeeze(distes(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end

% Ion analysis - OMNI
PSDiomni = zeros(length(disti.time),length(energyi0));
for ii = 1:length(disti.time);
    disttemp = squeeze(distis(ii,:,:,:));
    PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end

efluxomni = PSDomni.*energyspec.^2;
efluxomni = efluxomni/1e6/(5.486e-4)^2/0.53707; %convert to normal units

ifluxomni = PSDiomni.*energyspeci.^2;
ifluxomni = ifluxomni/1e6/0.53707; %convert to normal units

%% Compute PADS
[paddiste,thetae,energye,tinte] = mms.get_pitchangledist(diste,phie,thetae,stepTablee,energye0,energye1,Bxyz,tint);
[paddisti,thetai,energyi,tinti] = mms.get_pitchangledist(disti,phii,thetai,stepTablei,energyi0,energyi1,Bxyz,tint);

% Convert to DEflux
paddiste = paddiste/1e6/(5.486e-4)^2/0.53707;
paddisti = paddisti/1e6/0.53707;

for ii = 1:length(paddiste.time),
    energytemp = energye(ii,:)'*ones(1,length(thetae));
    paddiste.data(ii,:,:) = squeeze(paddiste.data(ii,:,:)).*energytemp.^2;
end

for ii = 1:length(paddisti.time),
    energytemp = energyi(ii,:)'*ones(1,length(thetae));
    paddisti.data(ii,:,:) = squeeze(paddisti.data(ii,:,:)).*energytemp.^2;
end

%%
E1 = [20 200];
E2 = [200 3000];
E3 = [3000 30000];

idxlow = find(energye0 > E1(1) & energye0 < E1(2));
idxmid = find(energye0 > E2(1) & energye0 < E2(2));
idxhigh= find(energye0 > E3(1) & energye0 < E3(2));

paddistelow = squeeze(mean(paddiste.data(:,idxlow,:),2));
paddistemid = squeeze(mean(paddiste.data(:,idxmid,:),2));
paddistehigh = squeeze(mean(paddiste.data(:,idxhigh,:),2));

paddistilow = squeeze(mean(paddisti.data(:,1:idxlow(end),:),2));
paddistimid = squeeze(mean(paddisti.data(:,idxmid,:),2));
paddistihigh = squeeze(mean(paddisti.data(:,idxhigh,:),2));

% Make structures for plotting
speceomni=struct('t',diste.time.epochUnix);
speceomni.p = double(efluxomni);
speceomni.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
speceomni.f_label={''};
speceomni.f = single(energyspec);

specepadl=struct('t',paddiste.time.epochUnix);
specepadl.p = double(paddistelow);
specepadl.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specepadl.f_label={''};
specepadl.f = single(thetae);

specepadm=struct('t',paddiste.time.epochUnix);
specepadm.p = double(paddistemid);
specepadm.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specepadm.f_label={''};
specepadm.f = single(thetae);

specepadh=struct('t',paddiste.time.epochUnix);
specepadh.p = double(paddistehigh);
specepadh.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specepadh.f_label={''};
specepadh.f = single(thetae);

speciomni=struct('t',disti.time.epochUnix);
speciomni.p = double(ifluxomni);
speciomni.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
speciomni.f_label={''};
speciomni.f = single(energyspeci);

specipadl=struct('t',paddisti.time.epochUnix);
specipadl.p = double(paddistilow);
specipadl.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
specipadl.f_label={''};
specipadl.f = single(thetai);

specipadm=struct('t',paddisti.time.epochUnix);
specipadm.p = double(paddistimid);
specipadm.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
specipadm.f_label={''};
specipadm.f = single(thetai);

specipadh=struct('t',paddisti.time.epochUnix);
specipadh.p = double(paddistihigh);
specipadh.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
specipadh.f_label={''};
specipadh.f = single(thetai);

%% Plot Ion data
h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)

h(2)=irf_panel('ni');
irf_plot(h(2),ni);
ylabel(h(2),'n_{i} (cm^{-3})','Interpreter','tex');
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)

h(3)=irf_panel('Vixyz');
irf_plot(h(3),Vi);
ylabel(h(3),'V_{i} (km s^{-1})','Interpreter','tex');
irf_legend(h(3),{'V_{x}','V_{y}','V_{z}'},[0.1 0.12])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)

h(4)=irf_panel('Ti');
irf_plot(h(4),Tiparperp);
ylabel(h(4),'T_{i} (eV)','Interpreter','tex');
irf_legend(h(4),{'T_{||}','T_{\perp}'},[0.1 0.12])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

h(5)=irf_panel('idist');
irf_spectrogram(h(5),speciomni,'log','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{i} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('ipadlow');
irf_spectrogram(h(6),specipadl,'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12);
set(h(6),'yscale','lin');
set(h(6),'ytick',[0 90 180]);
ylabel(h(6),{'\theta (deg.)','low E'},'fontsize',12,'Interpreter','tex');

h(7)=irf_panel('ipadmid');
irf_spectrogram(h(7),specipadm,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12);
set(h(7),'yscale','lin');
set(h(7),'ytick',[0 90 180]);
ylabel(h(7),{'\theta (deg.)','mid E'},'fontsize',12,'Interpreter','tex');

h(8)=irf_panel('ipadhigh');
irf_spectrogram(h(8),specipadh,'log','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12);
set(h(8),'yscale','lin');
set(h(8),'ytick',[0 90 180]);
ylabel(h(8),{'\theta (deg.)','high E'},'fontsize',12,'Interpreter','tex');

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);
colormap(h(7),cmap);
colormap(h(8),cmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',tinti);
set(h(1:8),'fontsize',12);

%% Plot electron data
h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)

h(2)=irf_panel('ne');
irf_plot(h(2),ne);
ylabel(h(2),'n_{e} (cm^{-3})','Interpreter','tex');
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)

h(3)=irf_panel('Vexyz');
irf_plot(h(3),Ve);
ylabel(h(3),'V_{e} (km s^{-1})','Interpreter','tex');
irf_legend(h(3),{'V_{x}','V_{y}','V_{z}'},[0.1 0.12])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)

h(4)=irf_panel('Te');
irf_plot(h(4),Teparperp);
ylabel(h(4),'T_{e} (eV)','Interpreter','tex');
irf_legend(h(4),{'T_{||}','T_{\perp}'},[0.1 0.12])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

h(5)=irf_panel('edist');
irf_spectrogram(h(5),speceomni,'log','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
hold(h(5),'on');
SCpot2 = SCpot.resample(Te);
Teav = TSeries(Te.time,[SCpot2.data (Te.data(:,1)+Te.data(:,4)+Te.data(:,6))/3]);
irf_plot(h(5),Teav)
hold(h(5),'off');
irf_legend(h(5),{'V_{SC}','T_{e}'},[0.38 0.98]);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('epadlow');
irf_spectrogram(h(6),specepadl,'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12);
set(h(6),'yscale','lin');
set(h(6),'ytick',[0 90 180]);
ylabel(h(6),{'\theta (deg.)','low E'},'fontsize',12,'Interpreter','tex');

h(7)=irf_panel('epadmid');
irf_spectrogram(h(7),specepadm,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12);
set(h(7),'yscale','lin');
set(h(7),'ytick',[0 90 180]);
ylabel(h(7),{'\theta (deg.)','mid E'},'fontsize',12,'Interpreter','tex');

h(8)=irf_panel('epadhigh');
irf_spectrogram(h(8),specepadh,'log','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12);
set(h(8),'yscale','lin');
set(h(8),'ytick',[0 90 180]);
ylabel(h(8),{'\theta (deg.)','high E'},'fontsize',12,'Interpreter','tex');

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);
colormap(h(7),cmap);
colormap(h(8),cmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',tinti);
set(h(1:8),'fontsize',12);