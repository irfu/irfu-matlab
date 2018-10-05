mmsId = 4;
Tint = irf.tint('2015-09-19T09:09:00Z/2015-09-19T09:11:00Z');

%% FGM & EDP
B_dmpa_fgm_srvy = mms.get_data('B_dmpa_fgm_srvy_l2',Tint,mmsId);
if isempty(B_dmpa_fgm_srvy)
  irf.log('warning','loading L2pre DFG')
  B_dmpa_fgm_srvy = mms.get_data('B_dmpa_dfg_srvy_l2pre',Tint,mmsId);
  if isempty(B_dmpa_fgm_srvy)
    irf.log('warning','loading QL DFG')
    B_dmpa_fgm_srvy = mms.get_data('B_dmpa_dfg_srvy_ql',Tint,mmsId);
  end
end
E_dsl_edp = mms.get_data('E_dsl_edp_brst_l2',Tint,mmsId);
if isempty(E_dsl_edp)
  irf.log('warning','loading QL DCE')
  E_dsl_edp = mms.get_data('E_dsl_edp_brst_ql',Tint,mmsId);
end
E2d_dsl_edp = mms.get_data('E2d_dsl_edp_brst_l2pre',Tint,mmsId);
if isempty(E2d_dsl_edp)
  irf.log('warning','loading QL DCE2d')
  E2d_dsl_edp = mms.get_data('E2d_dsl_edp_brst_ql',Tint,mmsId);
end
E_adp_edp = mms.get_data('E_ssc_edp_brst_l1b',Tint,mmsId);
E_adp_edp = -E_adp_edp.z*1.5;

%% FPI
Vi_dbcs_fpi = mms.get_data('Vi_dbcs_fpi_brst_l2',Tint,mmsId);
Ve_dbcs_fpi = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint,mmsId);
Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_brst_l2',Tint,mmsId);
if isempty(Vhplus_dbcs_hpca)
  Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_brst_l1b',Tint,mmsId);
end

%%
[~, Vi_perp] = irf_dec_parperp(B_dmpa_fgm_srvy,Vi_dbcs_fpi);
[~, Ve_perp] = irf_dec_parperp(B_dmpa_fgm_srvy,Ve_dbcs_fpi);
[~, E_perp] = irf_dec_parperp(B_dmpa_fgm_srvy,E_dsl_edp);
VExB = irf_e_vxb(E_dsl_edp,B_dmpa_fgm_srvy,-1);
VE2dxB = irf_e_vxb(E2d_dsl_edp,B_dmpa_fgm_srvy,-1);
EVixB = irf_e_vxb(Vi_dbcs_fpi,B_dmpa_fgm_srvy.resample(Vi_dbcs_fpi));
EVexB = irf_e_vxb(Ve_dbcs_fpi,B_dmpa_fgm_srvy.resample(Ve_dbcs_fpi));
if isempty(Vhplus_dbcs_hpca), Vhplus_perp = []; EVphlusxB = [];
else
  [~, Vhplus_perp] = irf_dec_parperp(B_dmpa_fgm_srvy,Vhplus_dbcs_hpca);
  EVphlusxB = irf_e_vxb(Vhplus_dbcs_hpca,B_dmpa_fgm_srvy.resample(Vhplus_dbcs_hpca));
end

%%
f = irf_figure(2387456,3);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.9 0;0 1 1 ;1 0 1; 1 1 0])
set(gcf,'defaultAxesFontSize',20)
h = irf_plot({VExB,VE2dxB,Ve_perp,Vi_perp,Vhplus_perp},'comp');
title(h(1),sprintf('MMS%d',mmsId))
legend(h(1),'VExB','VE2dxB','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
legend(h(2),'VExB','VE2dxB','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
legend(h(3),'VExB','VE2dxB','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
ylabel(h(1),'V_x DSL [km/s]')
ylabel(h(2),'V_y DSL [km/s]')
ylabel(h(3),'V_z DSL [km/s]')
irf_zoom(h,'x',Tint)

irf_print_fig(['mms' num2str(mmsId) '_VExB_EDP_vs_FPI_vs_HPCA_brst_' irf_fname(Tint,2)],'png')

%%
f = irf_figure(2387457,3);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.9 0;0 1 1 ;1 0 1; 1 1 0])
set(gcf,'defaultAxesFontSize',20)
h = irf_plot({E2d_dsl_edp,E_dsl_edp,E_perp,EVexB,EVixB,EVphlusxB},'comp');
irf_plot(h(3),E_adp_edp)
title(h(1),sprintf('MMS%d',mmsId))
legend(h(1),'E2d','E','E_perp','V_{e}xB','V_{i}xB','V_{H+}xB')
legend(h(2),'E2d','E','E_perp','V_{e}xB','V_{i}xB','V_{H+}xB')
legend(h(3),'E2d','E','E_perp','V_{e}xB','V_{i}xB','E adp')
ylabel(h(1),'E_x DSL [mV/m]')
ylabel(h(2),'E_y DSL [mV/m]')
ylabel(h(3),'E_z DSL [mV/m]')
irf_zoom(h,'x',Tint)

irf_print_fig(['mms' num2str(mmsId) '_E_EDP_vs_FPI_vs_HPCA_brst_' irf_fname(Tint,2)],'png')