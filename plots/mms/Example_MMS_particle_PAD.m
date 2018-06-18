% Calculate and plot electron and ion pitch angle distributions from L1b particle 
% brst data. 
% Written by D. B. Graham.

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

% Other variables
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_dmpa_srvy_l2'',Tint);',ic);
c_eval('Bgse=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',Tint);',ic);
c_eval('SCpot=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',Tint);',ic);

%% Compute parallel and perpendicular electron and ion temperatures
Tipp = mms.rotate_tensor(Ti,'fac',Bxyz,'pp'); 
Tiparperp = TSeries(Ti.time,[Tipp.xx.data Tipp.yy.data]);

Tepp = mms.rotate_tensor(Te,'fac',Bxyz,'pp'); 
Teparperp = TSeries(Te.time,[Tepp.xx.data Tepp.yy.data]);

%% Compute pitch-angle distributions and particle energy fluxes
% electron and ion omnidirection differential energy flux
ePDistomni = ePDist.omni.deflux; 
iPDistomni = iPDist.omni.deflux;

% electron and ion omnidirection differential energy flux
ePDistpitch = ePDist.e64.pitchangles(Bxyz,24).deflux;
iPDistpitch = iPDist.e64.pitchangles(Bxyz,18).deflux;

% select energy intervals for pitch-angle distributions
eint1 = [3e1 3e2];
eint2 = [3e2 3e3];
eint3 = [3e3 3e4];

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
irf_zoom(h(4),'y',[min(min(Tiparperp.data))/1.5 max(max(Tiparperp.data))*1.5]);
set(h(4),'yscale','log');
ylabel(h(4),'T_{i} (eV)','Interpreter','tex');
irf_legend(h(4),{'T_{||}','T_{\perp}'},[0.1 0.12])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

h(5)=irf_panel('idist');
irf_spectrogram(h(5),iPDistomni.specrec,'log','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{i} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('ipadlow');
irf_spectrogram(h(6),iPDistpitch.elim(eint1).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12);
set(h(6),'yscale','lin');
set(h(6),'ytick',[0 90 180]);
irf_legend(h(6),[num2str(eint1(1),'%g') 'eV -' num2str(eint1(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(6),{'\theta (deg.)','low E'},'fontsize',12,'Interpreter','tex');

h(7)=irf_panel('ipadmid');
irf_spectrogram(h(7),iPDistpitch.elim(eint2).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12);
set(h(7),'yscale','lin');
set(h(7),'ytick',[0 90 180]);
irf_legend(h(7),[num2str(eint2(1),'%g') 'eV -' num2str(eint2(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(7),{'\theta (deg.)','mid E'},'fontsize',12,'Interpreter','tex');

h(8)=irf_panel('ipadhigh');
irf_spectrogram(h(8),iPDistpitch.elim(eint3).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12);
set(h(8),'yscale','lin');
set(h(8),'ytick',[0 90 180]);
irf_legend(h(8),[num2str(eint3(1),'%g') 'eV -' num2str(eint3(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(8),{'\theta (deg.)','high E'},'fontsize',12,'Interpreter','tex');

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);
colormap(h(7),cmap);
colormap(h(8),cmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
set(h(1:8),'fontsize',12);

%% Plot electron data
h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bgse);
ylabel(h(1),'B (nT)','Interpreter','tex');
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
irf_zoom(h(4),'y',[min(min(Teparperp.data))/1.5 max(max(Teparperp.data))*1.5]);
set(h(4),'yscale','log');
ylabel(h(4),'T_{e} (eV)','Interpreter','tex');
irf_legend(h(4),{'T_{||}','T_{\perp}'},[0.1 0.12])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

h(5)=irf_panel('edist');
irf_spectrogram(h(5),ePDistomni.specrec,'log','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12);
hold(h(5),'on');
SCpot2 = SCpot.resample(Te);
Teav = TSeries(Te.time,[SCpot2.data Te.trace.data/3]);
irf_plot(h(5),Teav)
hold(h(5),'off');
irf_legend(h(5),{'V_{SC}','T_{e}'},[0.38 0.98]);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('epadlow');
irf_spectrogram(h(6),ePDistpitch.elim(eint1).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12);
set(h(6),'yscale','lin');
set(h(6),'ytick',[0 90 180]);
irf_legend(h(6),[num2str(eint1(1),'%g') 'eV -' num2str(eint1(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(6),{'\theta (deg.)','low E'},'fontsize',12,'Interpreter','tex');

h(7)=irf_panel('epadmid');
irf_spectrogram(h(7),ePDistpitch.elim(eint2).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12);
set(h(7),'yscale','lin');
set(h(7),'ytick',[0 90 180]);
irf_legend(h(7),[num2str(eint2(1),'%g') 'eV -' num2str(eint2(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(7),{'\theta (deg.)','mid E'},'fontsize',12,'Interpreter','tex');

h(8)=irf_panel('epadhigh');
irf_spectrogram(h(8),ePDistpitch.elim(eint3).specrec('pa'),'log','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12);
set(h(8),'yscale','lin');
set(h(8),'ytick',[0 90 180]);
irf_legend(h(8),[num2str(eint3(1),'%g') 'eV -' num2str(eint3(2),'%g') ' eV'],[0.02 0.1],'color','k')
ylabel(h(8),{'\theta (deg.)','high E'},'fontsize',12,'Interpreter','tex');

load('caa/cmap.mat');
colormap(h(5),cmap);
colormap(h(6),cmap);
colormap(h(7),cmap);
colormap(h(8),cmap);

c_eval('title(h(1),''MMS ?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
set(h(1:8),'fontsize',12);