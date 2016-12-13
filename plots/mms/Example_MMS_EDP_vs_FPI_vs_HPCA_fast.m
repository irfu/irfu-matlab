mmsId = 1;
Tint = irf.tint('2016-08-10T09:50:00Z/2016-08-10T10:15:00Z');

%% FGM & EDP
B_bcs_fgm_srvy_l2 = mms.get_data('B_dmpa_fgm_srvy_l2',Tint,mmsId);
E_dsl_edp_l2 = mms.get_data('E_dsl_edp_fast_l2',Tint,mmsId);
E2d_dsl_edp_l2pre = mms.get_data('E2d_dsl_edp_fast_l2pre',Tint,mmsId);

%% FPI
Vi_dbcs_fpi = mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,mmsId);
Ve_dbcs_fpi = mms.get_data('Ve_dbcs_fpi_fast_l2',Tint,mmsId);
Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_srvy_l2',Tint,mmsId);
if isempty(Vhplus_dbcs_hpca)
  Vhplus_dbcs_hpca = mms.get_data('Vhplus_dbcs_hpca_srvy_l1b',Tint,mmsId);
end

%%
[~, Vi_perp] = irf_dec_parperp(B_bcs_fgm_srvy_l2,Vi_dbcs_fpi);
[~, Ve_perp] = irf_dec_parperp(B_bcs_fgm_srvy_l2,Ve_dbcs_fpi);
VExB = irf_e_vxb(E_dsl_edp_l2,B_bcs_fgm_srvy_l2,-1);
VExB_l2pre = irf_e_vxb(E2d_dsl_edp_l2pre,B_bcs_fgm_srvy_l2,-1);
EVixB = irf_e_vxb(Vi_dbcs_fpi,B_bcs_fgm_srvy_l2.resample(Vi_dbcs_fpi));
EVexB = irf_e_vxb(Ve_dbcs_fpi,B_bcs_fgm_srvy_l2.resample(Ve_dbcs_fpi));
if isempty(Vhplus_dbcs_hpca), Vhplus_perp = []; EVphlusxB = [];
else
  [~, Vhplus_perp] = irf_dec_parperp(B_bcs_fgm_srvy_l2,Vhplus_dbcs_hpca);
  EVphlusxB = irf_e_vxb(Vhplus_dbcs_hpca,B_bcs_fgm_srvy_l2.resample(Vhplus_dbcs_hpca));
end

%%
f = irf_figure(2387456,3);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.7 0;0 1 1 ;1 0 1; 1 1 0])
h = irf_plot({VExB,VExB_l2pre,Ve_perp,Vi_perp,Vhplus_perp},'comp');
title(h(1),sprintf('MMS%d',mmsId))
legend(h(1),'VExB','VExB l2pre','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
legend(h(2),'VExB','VExB l2pre','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
legend(h(3),'VExB','VExB l2pre','V_{e\perp}','V_{i\perp}','V_{H+\perp}')
ylabel(h(1),'V_x DSL [km/s]')
ylabel(h(2),'V_y DSL [km/s]')
ylabel(h(3),'V_z DSL [km/s]')
irf_zoom(h,'x',Tint)

irf_print_fig(['mms' num2str(mmsId) '_VExB_EDP_vs_FPI_vs_HPCA_' irf_fname(Tint,2)],'png')

%%
f = irf_figure(2387457,3);
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.7 0;0 1 1 ;1 0 1; 1 1 0])
h = irf_plot({E2d_dsl_edp_l2pre,E_dsl_edp_l2,EVexB,EVixB,EVphlusxB},'comp');
title(h(1),sprintf('MMS%d',mmsId))
legend(h(1),'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB')
legend(h(2),'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB')
legend(h(3),'E L2pre','E l2','V_{e}xB','V_{i}xB','V_{H+}xB')
ylabel(h(1),'E_x DSL [mV/m]')
ylabel(h(2),'E_y DSL [mV/m]')
ylabel(h(3),'E_z DSL [mV/m]')
irf_zoom(h,'x',Tint)

irf_print_fig(['mms' num2str(mmsId) '_E_EDP_vs_FPI_vs_HPCA_' irf_fname(Tint,2)],'png')