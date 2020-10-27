%Example script to plot B, E, and ScPot (probe-to-spacecraft potential)

Tint = irf.tint('2020-06-24T12:00:00Z/2020-06-24T23:59:59Z');

%%
PSP = solo.db_get_ts('solo_L3_rpw-bia-scpot', 'PSP', Tint);
EDC_SRF = solo.db_get_ts('solo_L3_rpw-bia-efield', 'EDC_SRF', Tint);
EDC_RTN = EDC_SRF; EDC_RTN.data(:,1:2) = -EDC_RTN.data(:,1:2);
B_RTN = solo.db_get_ts('solo_L2_mag-rtn-normal','B_RTN', Tint);

%%
h = irf_figure(3094856,4,'reset');

hca = irf_panel('PSP');
irf_plot(hca,PSP);

hca = irf_panel('Btot');
irf_plot(hca,B_RTN.abs);
ylabel(hca,'|B| [nT]');

hca = irf_panel('Bxyz');
irf_plot(hca,B_RTN);
ylabel(hca,'B RTN [nT]');

hca = irf_panel('Exyz');
irf_plot(hca,EDC_SRF);
hca.YLim = [-9 9];
ylabel(hca,'E RTN [mV/m]');

irf_zoom(h,'x',Tint)