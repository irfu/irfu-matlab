mmsId = 1;
Tint = irf.tint('2016-08-10T09:50:00Z/2016-08-10T10:15:00Z');
fpiMode = 'fast'; % alternative fpiMode = 'brst'
edpMode = 'fast'; % alternative edpMode = 'brst'

%% FGM & EDP
switch edpMode
  case 'fast', fgmMode = 'srvy';
  case 'brst', fgmMode = 'brst';
  otherwise, error('unrecognized mode (FAST/BRST)')
end
B_dmpa_fgm_srvy_l2 = mms.get_data(['B_dmpa_fgm_' fgmMode '_l2'],Tint,mmsId);
E_dsl_edp_l2 = mms.get_data(['E_dsl_edp_' edpMode '_l2'],Tint,mmsId);
E2d_dsl_edp_l2pre = mms.get_data(['E2d_dsl_edp_' edpMode '_l2pre'],Tint,mmsId);
E_adp_edp = mms.get_data(['E_ssc_edp_' edpMode '_l1b'],Tint,mmsId);
E_adp_edp = -E_adp_edp.z*1.5;

%% FPI
fpiSuf = ['_fpi_' fpiMode '_l2'];
Vi_dbcs_fpi = mms.get_data(['Vi_dbcs' fpiSuf],Tint,mmsId);
Ve_dbcs_fpi = mms.get_data(['Ve_dbcs' fpiSuf],Tint,mmsId);
Ne_fpi = mms.get_data(['Ne' fpiSuf],Tint,mmsId);
idx_lowDensity = Ne_fpi.data < 0.06; Ve_dbcs_fpi.data(idx_lowDensity,:) = NaN;
Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_srvy_l2',Tint,mmsId);
if isempty(Vhplus_dbcs_hpca)
  Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_srvy_l1b',Tint,mmsId);
end

% correct Ez in E2d
% XXX: this should be undone later
E2d_dsl_edp_l2pre = irf_edb(E2d_dsl_edp_l2pre,B_dmpa_fgm_srvy_l2,10,'Eperp+NaN');

% Comp VxB
[~, Vi_perp] = irf_dec_parperp(B_dmpa_fgm_srvy_l2,Vi_dbcs_fpi);
[~, Ve_perp] = irf_dec_parperp(B_dmpa_fgm_srvy_l2,Ve_dbcs_fpi);
VExB = irf_e_vxb(E_dsl_edp_l2,B_dmpa_fgm_srvy_l2,-1);
VExB_l2pre = irf_e_vxb(E2d_dsl_edp_l2pre,B_dmpa_fgm_srvy_l2,-1);
EVixB = irf_e_vxb(Vi_dbcs_fpi,B_dmpa_fgm_srvy_l2.resample(Vi_dbcs_fpi));
EVexB = irf_e_vxb(Ve_dbcs_fpi,B_dmpa_fgm_srvy_l2.resample(Ve_dbcs_fpi));
if isempty(Vhplus_dbcs_hpca), Vhplus_perp = []; EVphlusxB = [];
else
  [~, Vhplus_perp] = irf_dec_parperp(B_dmpa_fgm_srvy_l2,Vhplus_dbcs_hpca);
  EVphlusxB = irf_e_vxb(Vhplus_dbcs_hpca,B_dmpa_fgm_srvy_l2.resample(Vhplus_dbcs_hpca));
end

%% E plot
h = irf_figure(2388456,4);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.7 0;0 1 1 ;1 0 1; 1 1 0])
set(gcf,'defaultAxesFontSize',12)

hca = irf_panel('b');
irf_plot(hca,B_dmpa_fgm_srvy_l2)
title(hca,sprintf('MMS%d',mmsId))
ylabel(hca,'B DSL [nT]')
irf_legend(hca,{'X','Y','Z'},[0.98 0.05],'fontsize',14)

hca = irf_panel('Ex');
irf_plot(hca,{E2d_dsl_edp_l2pre.x,E_dsl_edp_l2.x,EVexB.x,EVixB.x,EVphlusxB.x},'comp');
irf_legend(hca,{'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'E_x DSL [mV/m]')

hca = irf_panel('Ey');
irf_plot(hca,{E2d_dsl_edp_l2pre.y,E_dsl_edp_l2.y,EVexB.y,EVixB.y,EVphlusxB.y},'comp');
irf_legend(hca,{'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'E_y DSL [mV/m]')

hca = irf_panel('Ez');
irf_plot(hca,...
  {E2d_dsl_edp_l2pre.z,E_dsl_edp_l2.z,EVexB.z,EVixB.z,EVphlusxB.z,E_adp_edp},'comp');
irf_legend(hca,{'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB','E ADP'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'E_z DSL [mV/m]')

irf_zoom(h,'x',Tint)
irf_plot_ylabels_align(h)
irf_plot_axis_align(h)

irf_print_fig(['mms' num2str(mmsId) '_E_EDP_' edpMode '_vs_FPI_' fpiMode '_vs_HPCA_fast_' irf_fname(Tint,2)],'png')

%% V plot
h = irf_figure(2387458,4);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.7 0;0 1 1 ;1 0 1; 1 1 0])
set(gcf,'defaultAxesFontSize',12)

hca = irf_panel('b');
irf_plot(hca,B_dmpa_fgm_srvy_l2)
title(hca,sprintf('MMS%d',mmsId))
ylabel(hca,'B DSL [nT]')
irf_legend(hca,{'X','Y','Z'},[0.98 0.05],'fontsize',14)

hca = irf_panel('Vx');
irf_plot(hca,{VExB_l2pre.x,VExB.x,Ve_perp.x,Vi_perp.x,Vhplus_perp.x},'comp');
irf_legend(hca,{'VExB l2pre','VExB','V_{e\perp}','V_{i\perp}','V_{H+\perp}'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'V_x DSL [km/s]')

hca = irf_panel('Vy');
irf_plot(hca,{VExB_l2pre.y,VExB.y,Ve_perp.y,Vi_perp.y,Vhplus_perp.y},'comp');
irf_legend(hca,{'VExB l2pre','VExB','V_{e\perp}','V_{i\perp}','V_{H+\perp}'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'V_y DSL [km/s]')

hca = irf_panel('Vz');
irf_plot(hca,{VExB_l2pre.z,VExB.z,Ve_perp.z,Vi_perp.z,Vhplus_perp.z},'comp');
irf_legend(hca,{'VExB l2pre','VExB','V_{e\perp}','V_{i\perp}','V_{H+\perp}'},...
  [0.98 0.05],'fontsize',14)
ylabel(hca,'V_z DSL [km/s]')

irf_zoom(h,'x',Tint)
irf_plot_ylabels_align(h)
irf_print_fig(['mms' num2str(mmsId) '_VExB_EDP_' edpMode '_vs_FPI_' fpiMode '_vs_HPCA_fast_' irf_fname(Tint,2)],'png')


%%
if 0
  f = irf_figure(2387458,3);
  set(gcf,'defaultAxesColorOrder',[1 0 0;0 0 0;0 0 1;0 0.7 0;0 1 1 ;1 0 1; 1 1 0])
  set(gcf,'defaultAxesFontSize',12)
  h = irf_plot({EVexB,E2d_dsl_edp_l2pre,E_dsl_edp_l2,EVixB,EVphlusxB},'comp');
  title(h(1),sprintf('MMS%d',mmsId))
  legend(h(1),'V_{e}xB','E L2pre','E l2','V_{i}xB','V_{H+}xB','Location','NorthEastOutside')
  legend(h(2),'V_{e}xB','E L2pre','E l2','V_{i}xB','V_{H+}xB','Location','NorthEastOutside')
  legend(h(3),'V_{e}xB','E L2pre','E l2','V_{i}xB','V_{H+}xB','Location','NorthEastOutside')
  ylabel(h(1),'E_x DSL [mV/m]')
  ylabel(h(2),'E_y DSL [mV/m]')
  ylabel(h(3),'E_z DSL [mV/m]')
  irf_zoom(h,'x',Tint)
  irf_plot_ylabels_align(h)

  irf_print_fig(['mms' num2str(mmsId) '_E_EDP_vs_FPI_vs_HPCA_fast_' irf_fname(Tint,2)],'png')
end
