% Load brst L1b particle distributions and convert to differential energy
% fluxes. Plots electron and ion fluxes and electron anisotropies. 
%
% Written by D. B. Graham.
%

ic = 3; % Spacecraft number

Tint = irf.tint('2015-10-30T05:15:20.00Z/2015-10-30T05:16:20.00Z');

%% Load Data
% Particle distributions
c_eval('ePDist = mms.get_data(''PDe_fpi_brst_l2'',Tint,?);',ic)
c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',ic)

% Particle moments
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);
c_eval('Ve = mms.get_data(''Ve_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('Te = mms.get_data(''Te_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('ni = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',Tint);',ic);
c_eval('Vi = mms.get_data(''Vi_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('Ti = mms.get_data(''Ti_gse_fpi_brst_l2'',Tint,?);',ic);

c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',Tint);',ic);
c_eval('Bgse=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);

c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
c_eval('SCpot=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',Tint);',ic);

%% Compute moments
Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
mp = Units.mp;

Vav = Vi.abs.data*1000;
Vav = 0.5*mp*Vav.^2/qe;
Vav = irf.ts_scalar(Vi.time,Vav);

ePDistomni = ePDist.omni.deflux; 
iPDistomni = iPDist.omni.deflux;

ePDistpitch = ePDist.pitchangles(Bxyz,13).deflux;

% Make structures for plotting
PSDparapar = ePDistpitch.data(:,:,1)./ePDistpitch.data(:,:,13);
specparapar =struct('t',ePDistpitch.time.epochUnix);
specparapar.p = PSDparapar;
specparapar.p_label={'log_{10}(f_{||+}/f_{||-})'};
specparapar.f_label={''};
specparapar.f = ePDistomni.depend{1,1};

PSDparperp = (ePDistpitch.data(:,:,1)+ePDistpitch.data(:,:,13))./(2*ePDistpitch.data(:,:,7));
specparperp =struct('t',ePDistpitch.time.epochUnix);
specparperp.p = PSDparperp;
specparperp.p_label={'log_{10}(f_{||+}+f_{||-})/(2 f_{\perp})'};
specparperp.f_label={''};
specparperp.f = ePDistomni.depend{1,1};

%% Plot Figure

h=irf_plot(8,'newfigure');
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bgse);
ylabel(h(1),'B (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)

h(2)=irf_panel('Vxyz');
irf_plot(h(2),Vi);
ylabel(h(2),'V_{i} (km s^{-1})','Interpreter','tex');
irf_legend(h(2),{'V_{x}','V_{y}','V_{z}'},[0.1 0.12])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)

ni = ni.resample(ne);
nall = TSeries(ne.time,[ne.data ni.data]);

h(3)=irf_panel('ne');
irf_plot(h(3),nall);
set(h(3),'yscale','log');
set(h(3),'ytick',[1e-1 1e0 1e1 1e2]);
irf_zoom(h(3),'y',[min(min(nall.data))/1.5 max(max(nall.data))*1.5]);
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
irf_spectrogram(h(5),iPDistomni.specrec,'log','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{i} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('edist');
irf_spectrogram(h(6),ePDistomni.specrec,'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12);
hold(h(6),'on');
SCpot2 = SCpot.resample(Te);
Teav = TSeries(Te.time,[SCpot2.data Te.trace.data/3]);
irf_plot(h(6),Teav)
hold(h(6),'off');
irf_legend(h(6),{'V_{SC}','T_{e}'},[0.38 0.98]);
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

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);
colormap(h(7),bgrcmap);
colormap(h(8),bgrcmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
set(h(1:8),'fontsize',12);