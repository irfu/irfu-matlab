%% Set time
Tint = irf.tint('2020-08-28T06:00:00Z/2020-08-28T12:59:59Z');

%% Load data
B_SRF = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);
V_SRF = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);
Np = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);
Ne = solo.db_get_ts('solo_L3_rpw-bia-density', 'DENSITY', Tint);
EDC = solo.db_get_ts('solo_L3_rpw-bia-efield', 'EDC_SRF', Tint);

%% Plot
h=irf_plot(5,'newfigure');
%
hca=irf_panel('B');
irf_plot(hca,B_SRF.abs);
ylabel(hca,'B [nT]');
%
hca=irf_panel('BXYZ');
irf_plot(hca,B_SRF);
ylabel(hca,'B [nT] SRF');
irf_legend(hca,{'B_X';'B_Y';'B_Z'},[1.02 0.98])
%
hca=irf_panel('n');
irf_plot(hca,{Ne,Np},'comp');
ylabel(hca,'n [cc]');
irf_legend(hca,{'Ne';'Np'},[1.02 0.98])
%
hca=irf_panel('V');
irf_plot(hca,V_SRF);
ylabel(hca,'V [km/s] SRF');
irf_legend(hca,{'V_X';'V_Y';'V_Z'},[1.02 0.98])
%
hca=irf_panel('E');
irf_plot(hca,EDC);
ylabel(hca,'E [mV/m] SRF');
irf_legend(hca,{'E_X';'E_Y';'E_Z'},[1.02 0.98])
%
axis(h,'tight')
%set(hca,'YLim',19*[-1 1])
irf_plot_axis_align(h)