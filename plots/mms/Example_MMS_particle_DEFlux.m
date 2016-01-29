% Load brst L1b particle distributions and convert to differential energy
% fluxes. Plots electron and ion fluxes and electron anisotropies. 
%
% Written by D. B. Graham.
%

ic = 1; % Spacecraft number

tint = irf.tint('2015-12-30T00:30:00.00Z/2015-12-30T00:30:35.00Z');
%%
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

c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
Bxyz = Bxyz.resample(diste);
Bvec = Bxyz/Bxyz.abs;

c_eval('Exyz=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

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
Vav = Vi.abs.data*1000;
Vav = 0.5*mp*Vav.^2/qe;
Vav = irf.ts_scalar(Vi.time,Vav);
emoments = mms.psd_moments(diste,phie,thetae,stepTablee,energye0,energye1,SCpot,'electron');
ne = emoments.n_psd;
Ve = emoments.V_psd;
Te = emoments.T_psd;

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

xe = zeros(length(diste.time),lengthphi,length(thetae));
ye = zeros(length(diste.time),lengthphi,length(thetae));
ze = -ones(lengthphi,1)*cosd(thetae);

for ii = 1:length(disti.time);
    for jj=1:length(energyi0);
        allsolidi(ii,jj,:,:) = solida;
    end
end

for ii = 1:length(diste.time);
    for jj=1:length(energye0);
        allsolide(ii,jj,:,:) = solida;
    end
    xe(ii,:,:) = -cosd(phie.data(ii,:)')*sind(thetae);
    ye(ii,:,:) = -sind(phie.data(ii,:)')*sind(thetae);
end

distis = disti.data.*allsolidi;
distes = diste.data.*allsolide;

% Define new PAD arrays
PSDomni = zeros(length(diste.time),length(energye0));
PSDpar = PSDomni;
PSDperp = PSDomni;
PSDapar = PSDomni;
PSDpartemp = zeros(length(energye0),lengthphi,length(thetae));
PSDperptemp = zeros(length(energye0),lengthphi,length(thetae));
PSDapartemp = zeros(length(energye0),lengthphi,length(thetae));

% Electron analysis
for ii = 1:length(diste.time);
    disttemp = squeeze(distes(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
    disttemp = squeeze(diste.data(ii,:,:,:));
    Bvecs = Bvec.data(ii,:);
    x = squeeze(xe(ii,:,:));
    y = squeeze(ye(ii,:,:));
    thetab = acosd(x*Bvecs(1)+y*Bvecs(2)+ze*Bvecs(3));
    pospar = ones(lengthphi,length(thetae)); 
    pospar(thetab > 15) = NaN;
    posperp = ones(lengthphi,length(thetae)); 
    posperp(thetab < 82.5) = NaN;
    posperp(thetab > 97.5) = NaN;
    posapar = ones(lengthphi,length(thetae)); 
    posapar(thetab < 165) = NaN; 
    for kk = 1:length(energye0);
        PSDpartemp(kk,:,:)  = squeeze(disttemp(kk,:,:)).*pospar;
        PSDperptemp(kk,:,:) = squeeze(disttemp(kk,:,:)).*posperp;
        PSDapartemp(kk,:,:) = squeeze(disttemp(kk,:,:)).*posapar;
    end
    PSDpar(ii,:) =  squeeze(irf.nanmean(irf.nanmean(PSDpartemp,3),2));
    PSDperp(ii,:) = squeeze(irf.nanmean(irf.nanmean(PSDperptemp,3),2));
    PSDapar(ii,:) = squeeze(irf.nanmean(irf.nanmean(PSDapartemp,3),2));
end

% Ion analysis
PSDiomni = zeros(length(disti.time),length(energyi0));
for ii = 1:length(disti.time);
    disttemp = squeeze(distis(ii,:,:,:));
    PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end

efluxomni = PSDomni.*energyspec.^2;
efluxomni = efluxomni/1e6/(5.486e-4)^2/0.53707; %convert to normal units

ifluxomni = PSDiomni.*energyspeci.^2;
ifluxomni = ifluxomni/1e6/0.53707; %convert to normal units

specomni=struct('t',diste.time.epochUnix);
specomni.p = double(PSDomni)*1e30;
specomni.p_label={'f_e','(s^{3} km^{-6})'};
specomni.f_label={''};
specomni.f = single(energyspec);

specfomni=struct('t',diste.time.epochUnix);
specfomni.p = double(efluxomni);
specfomni.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specfomni.f_label={''};
specfomni.f = single(energyspec);

specfiomni=struct('t',disti.time.epochUnix);
specfiomni.p = double(ifluxomni);
specfiomni.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
specfiomni.f_label={''};
specfiomni.f = single(energyspeci);

specpar =struct('t',diste.time.epochUnix);
specpar.p = double(PSDpar);
specpar.p_label={'f_e par','(s^{3} m^{-6})'};
specpar.f_label={''};
specpar.f = single(energyspec);

specapar =struct('t',diste.time.epochUnix);
specapar.p = double(PSDapar);
specapar.p_label={'f_e perp','(s^{3} m^{-6})'};
specapar.f = single(energyspec);

specperp =struct('t',diste.time.epochUnix);
specperp.p = double(PSDperp);
specperp.p_label={'f_e apar','(s^{3} m^{-6})'};
specperp.f_label={''};
specperp.f = single(energyspec);

PSDparapar = PSDpar./PSDapar;
specparapar =struct('t',diste.time.epochUnix);
specparapar.p = double(PSDparapar);
specparapar.p_label={'f_{par}/f_{apar}'};
specparapar.f_label={''};
specparapar.f = single(energyspec);

PSDparperp = (PSDpar+PSDapar)./(2*PSDperp);
specparperp =struct('t',diste.time.epochUnix);
specparperp.p = double(PSDparperp);
specparperp.p_label={'(f_{par}+f_{apar})/(2 f_{perp})'};
specparperp.f_label={''};
specparperp.f = single(energyspec);

%% Plot Figure

h=irf_plot(8,'newfigure');
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)

h(2)=irf_panel('Vxyz');
irf_plot(h(2),Vi);
ylabel(h(2),'V_{DSL} (km s^{-1})','Interpreter','tex');
irf_legend(h(2),{'V_{x}','V_{y}','V_{z}'},[0.1 0.12])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)

ni = ni.resample(ne);
nall = TSeries(ne.time,[ne.data ni.data]);

h(3)=irf_panel('ne');
irf_plot(h(3),nall);
set(h(3),'yscale','log');
set(h(3),'ytick',[1e-1 1e0 1e1 1e2]);
irf_zoom(h(3),'y',[0.1 50]);
irf_legend(h(3),{'n_{e}','n_{i}'},[0.10 0.12])
ylabel(h(3),'n (cm^{-3})','Interpreter','tex');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)

SCpot = SCpot.resample(ne);

h(4) = irf_panel('Exyz');
irf_plot(h(4),Exyz);
ylabel(h(4),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(4),{'E_{x}','E_{y}','E_{z}'},[0.10 0.12])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

h(5)=irf_panel('idist');
irf_spectrogram(h(5),specfiomni,'log');
hold(h(5),'on');
irf_plot(h(5),Vav);
hold(h(5),'off');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
irf_legend(h(5),{'V'},[0.38 0.4]);
ylabel(h(5),'E_{i} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('edist');
irf_spectrogram(h(6),specfomni,'log');
hold(h(6),'on');
Teav = TSeries(Te.time,[SCpot.data (Te.data(:,1)+Te.data(:,4)+Te.data(:,6))/3]);
irf_plot(h(6),Teav)
hold(h(6),'off');
irf_legend(h(6),{'V_{SC}','T_{e}'},[0.38 0.98]);
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(7)=irf_panel('edistparapar');
irf_spectrogram(h(7),specparapar,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)
hold(h(7),'on');
irf_plot(h(7),SCpot);
hold(h(7),'off');
set(h(7),'yscale','log');
caxis(h(7),[-2 2])
set(h(7),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(7),'E_{e} (eV)','fontsize',12);

h(8)=irf_panel('edistparperp');
irf_spectrogram(h(8),specparperp,'log','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12)
hold(h(8),'on');
irf_plot(h(8),SCpot);
hold(h(8),'off');
set(h(8),'yscale','log');
caxis(h(8),[-2 2])
set(h(8),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(8),'E_{e} (eV)','fontsize',12);

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);

load('caa/rgbcmap.mat');
colormap(h(7),rgbcmap);
colormap(h(8),rgbcmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',tint);
set(h(1:8),'fontsize',12);

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','MMS1overviewrev1.png');