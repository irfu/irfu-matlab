%%  remove electron distribution function background.
%

%%  1. parameter
clear;
ic = 1;
ie = 'e';
%    Tint = irf.tint('2016-12-24T15:03:10.000Z/2016-12-24T15:03:34.000Z');
%    time = irf_time('2016-12-24T15:03:32.054Z','utc>epochtt');
Tint = irf.tint('2015-12-02T01:14:40.000Z/2015-12-02T01:15:05.000Z');
time = irf_time('2015-12-02T01:15:00.000Z','utc>epochtt');
Npa = 24;
dir1 = [1, 0, 0];
dir2 = [0, 1, 0];
dir3 = [0, 0, 1];
dirstr = {'V_X', 'V_Y', 'V_Z'};

%%  2. load data
% 2.1. FPI
c_eval('Ne? = mms.get_data(''Ne_fpi_brst_l2'',Tint, ?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',Tint, ?);',ic);
c_eval('eDist? = mms.get_data(''PDe_fpi_brst_l2'',Tint, ?);',ic);
% 2.2. field
c_eval('gseB? = mms.get_data(''B_gse_fgm_brst_l2'', Tint, ?);',ic);
c_eval('dmpaB? = mms.get_data(''B_dmpa_fgm_brst_l2'', Tint, ?);',ic);
c_eval('gseE? = mms.get_data(''E_gse_edp_brst_l2'', Tint, ?);',ic);
c_eval('dslE? = mms.get_data(''E_dsl_edp_brst_l2'', Tint, ?);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',Tint);',ic);

%%  3. remove data for Tint
[eDist1_bgrm, eDist1_bg, nphotoe]= remove_edist_background(eDist1, 'tint', Tint, 'ZeroNaN', 0);

%%  4. Partial moments
if 0
  tic; %#ok<UNRCH>
  c_eval('epdist = eDist?;', ic);
  edistin = TSeries(epdist.time, epdist.data);
  ephi = TSeries(epdist.time, epdist.depend{1,2}(:,:));
  etheta = epdist.depend{1,3};
  estepTable = TSeries(epdist.time, epdist.ancillary.esteptable(:));
  eenergy0 = epdist.ancillary.energy0;
  eenergy1 = epdist.ancillary.energy1;
  energy_plus_tmp = epdist.ancillary.delta_energy_plus(1, :);
  energy_minus_tmp = epdist.ancillary.delta_energy_minus(1, :);
  epsd2mom_Tmp1 = mms.psd_moments(edistin, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32]);
  epsd2mom_Tmp2 = mms.psd_moments(edistin, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp, 'innerelec', 'on');
  epsd2mom_Tmp3 = mms.psd_moments(edistin, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp);
  epsd2mom_Tmp4 = mms.psd_moments(edistin, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [2, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp);
  %
  epsd2mom_tmp1 = mms.psd_moments(eDist1_bgrm, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32]);
  epsd2mom_tmp2 = mms.psd_moments(eDist1_bgrm, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp, 'innerelec', 'on');
  epsd2mom_tmp3 = mms.psd_moments(eDist1_bgrm, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [1, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp);
  epsd2mom_tmp4 = mms.psd_moments(eDist1_bgrm, ephi, etheta, estepTable, eenergy0, eenergy1, scPot1, 'e', 'enchannels', [2, 32], 'energy_plus', energy_plus_tmp, 'energy_minus', energy_minus_tmp);
  %
  nskymap = epsd2mom_Tmp3.n_psd_skymap;
  c_eval('nskymap.species = eDist?.species;', ic);
  c_eval('nskymap.ancillary.delta_energy_minus = eDist?.ancillary.delta_energy_minus;', ic);
  c_eval('nskymap.ancillary.delta_energy_plus = eDist?.ancillary.delta_energy_plus;', ic);
  %c_eval('nskymap.ancillary.energy = eDist?.ancillary.energy(1, :);', ic);
  c_eval('nskymap.ancillary.energy0 = eDist?.ancillary.energy0;', ic);
  c_eval('nskymap.ancillary.energy1 = eDist?.ancillary.energy1;', ic);
  c_eval('nskymap.ancillary.esteptable = eDist?.ancillary.esteptable(:, 1);', ic);
  ne32 = epsd2mom_Tmp3.n_psd_e32;
  toc;
end

%%  5. pitch angle distribution
Units = irf_units;
qe = Units.e;           me = Units.me;
c_eval('eDist_pad = eDist?.pitchangles(dmpaB?, Npa);', ic);
c_eval('eDist_bgrm_pad = eDist?_bgrm.pitchangles(dmpaB?, Npa);', ic);
c_eval('eDist_bg_pad = eDist?_bg.pitchangles(dmpaB?, Npa);', ic);
%%
c_eval('tdiff = mean(diff(eDist?.time.epochUnix));', ic);
c_eval('[~, nid] = min(abs(eDist?.time.epochUnix - time.epochUnix));', ic);
c_eval('time0 = eDist?.time(nid);', ic);
tint = irf.tint(time0 + (-tdiff/2.), time0 + (tdiff/2.));
c_eval('scPot_tint = scPot?.tlim(tint);', ic);
c_eval('scPot_tint = irf.nanmean(scPot_tint.data, 1);', ic);
c_eval('gseB?tint = gseB?.tlim(tint);', ic);
c_eval('gseB?tint = irf.nanmean(gseB?tint.data, 1);', ic);
c_eval('hatB? = gseB?tint/norm(gseB?tint);', ic);
% eDist @ tint
eDist_pad_tmp = eDist_pad.tlim(tint);
eDist_pad_tint = squeeze(irf.nanmean(eDist_pad_tmp.data, 1));
eDist_pad_pa = eDist_pad_tmp.depend{2};
energy_tint = eDist_pad_tmp.depend{1};
Energy_tint = energy_tint - scPot_tint;
Venergy_tint = sqrt(Energy_tint * qe * 2 / me);
eDist_para_tmp = irf.nanmean(eDist_pad_tint(:, 1:2)*1e30, 2);
eDist_perp_tmp = irf.nanmean(eDist_pad_tint(:, 12:13)*1e30, 2);
eDist_apar_tmp = irf.nanmean(eDist_pad_tint(:, 23:24)*1e30, 2);
% eDist_bgrm @ tint
eDist_bgrm_pad_tmp = eDist_bgrm_pad.tlim(tint);
eDist_bgrm_pad_tint = squeeze(irf.nanmean(eDist_bgrm_pad_tmp.data, 1));
eDist_bgrm_para_tmp = irf.nanmean(eDist_bgrm_pad_tint(:, 1:2)*1e30, 2);
eDist_bgrm_perp_tmp = irf.nanmean(eDist_bgrm_pad_tint(:, 12:13)*1e30, 2);
eDist_bgrm_apar_tmp = irf.nanmean(eDist_bgrm_pad_tint(:, 23:24)*1e30, 2);
% eDist_bg @ tint
eDist_bg_pad_tmp = eDist_bg_pad.tlim(tint);
eDist_bg_pad_tint = squeeze(irf.nanmean(eDist_bg_pad_tmp.data, 1));
eDist_bg_para_tmp = irf.nanmean(eDist_bg_pad_tint(:, 1:2)*1e30, 2);
eDist_bg_perp_tmp = irf.nanmean(eDist_bg_pad_tint(:, 12:13)*1e30, 2);
eDist_bg_apar_tmp = irf.nanmean(eDist_bg_pad_tint(:, 23:24)*1e30, 2);

%%  4. plot
% 4.0. figure setting
npanel1 = 7;    npanel2 = 3;       npanel3 = 1;
h = irf_plot(npanel1 + npanel2 + npanel3, 'reset');
xSize=1270; ySize=650;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth1 = 0.37;         ywidth1 = 0.9/npanel1;     xpos1 = 0.06;           ypos1 = 0.97;
xwidth2 = 0.25;         ywidth2 = 0.25;     xpos2 = 0.42;           ypos2 = 0.73;
xwidth3 = 0.40;         ywidth3 = 0.54;     xpos3 = 0.52;           ypos3 = 0.09;
c_eval('set(h(?), ''position'', [0.07 0.970-? * ywidth1 xwidth1 ywidth1]);', 1 : npanel1);

% 4.1. Bxyz
ip = 1;
h(ip) = irf_panel('Bgse');
irf_plot(h(ip), gseB1, 'LineWidth', 1.5);
%hold(h(ip), 'on');
%irf_plot(h(ip), lmnB.abs, 'LineWidth', 1.5);
%hold(h(ip), 'off');
ylabel(h(ip),{'B','[nT]'},'Interpreter','tex');
%irf_legend(h(ip),{['B_{' xyz(1) '} '],[' B_{' xyz(2) '} '], [' B_{' xyz(3) '} '], ' |B|'},[0.88 0.15])
irf_legend(h(ip),'(a)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');

% 5.L2. Ne
ip = ip + 1;
h(ip) = irf_panel('ne');
c_eval('irf_plot(h(ip), Ne?, ''LineWidth'', 1.5);', ic);
ylabel(h(ip),{'N_e', '[cm^{-3}]'},'Interpreter','tex');
irf_legend(h(ip),{'(b)'},[0.98 0.98], 'fontsize',12, 'fontWeight', 'bold' )


% 5.L3. gseVe
ip = ip + 1;
h(ip) = irf_panel('gseVe');
c_eval('irf_plot(h(ip), gseVe?, ''LineWidth'', 1.5);', ic);
ylabel(h(ip),{'V_e','[km/s]'},'Interpreter','tex');
irf_legend(h(ip),{'(c)'},[0.98 0.98], 'fontsize',12, 'fontWeight', 'bold' )
irf_legend(h(ip),{'V_{X} ',' V_{Y} ',' V_{Z}'},[0.98 0.12], 'fontsize',14, 'fontWeight', 'bold' )


% 4. gseE
ip = ip + 1;
h(ip) = irf_panel('gseE');
c_eval('irf_plot(h(ip), gseE?, ''LineWidth'', 1.5);', ic);
%hold(h(ip), 'on');
%c_eval('irf_plot(h(ip), gseVexB?.tlim(tint_tmp1).x, ''ko'', ''LineWidth'', 1.5, ''MarkerFacecolor'', ''k'');', ic);
%c_eval('irf_plot(h(ip), gseVexB?.tlim(tint_tmp1).y, ''bo'', ''LineWidth'', 1.5, ''MarkerFacecolor'', ''b'');', ic);
%c_eval('irf_plot(h(ip), gseVexB?.tlim(tint_tmp1).z, ''ro'', ''LineWidth'', 1.5, ''MarkerFacecolor'', ''r'');', ic);
%hold(h(ip), 'off');
ylabel(h(ip),{'E','[mV/m]'},'Interpreter','tex');
irf_legend(h(ip),{'(e)'},[0.98 0.98], 'color', 'k', 'fontsize',12, 'fontWeight', 'bold' )
irf_legend(h(ip),{'E_{X} ',' E_{Y}',' E_{Z}'},[0.02 0.15], 'fontsize',14, 'fontWeight', 'bold' )

% 4.7. electron energy flux
ip = ip + 1;
h(ip) = irf_panel('eEnflux_original');
c_eval('irf_spectrogram(h(ip), eDist?.deflux.omni.specrec(''energy''), ''log'');', ic);
hold(h(ip), 'on');
c_eval('irf_plot(h(ip), scPot?, ''LineWidth'', 1.3);', ic);
hold(h(ip), 'off');
h(ip).YScale = 'log';
h(ip).YTick = 10.^[1 2 3 4];
irf_legend(h(ip),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
%irf_legend(h(ip),'T_{||}',[0.75 0.1],'color','k')
%irf_legend(h(ip),'T_{\perp}',[0.8 0.1],'color','w')
ylabel(h(ip),{'W_e','[eV]'},'Interpreter','tex');
grid(h(ip), 'off');
caxis(h(ip), [3.5 9]);

% 4.7. electron energy flux
ip = ip + 1;
h(ip) = irf_panel('eEnflux_new');
c_eval('irf_spectrogram(h(ip), eDist?_bgrm.deflux.omni.specrec(''energy''), ''log'');', ic);
hold(h(ip), 'on');
c_eval('irf_plot(h(ip), scPot?, ''LineWidth'', 1.3);', ic);
hold(h(ip), 'off');
h(ip).YScale = 'log';
h(ip).YTick = 10.^[1 2 3 4];
caxis(h(ip), [3.5 9]);
irf_legend(h(ip),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
%irf_legend(h(ip),'T_{||}',[0.75 0.1],'color','k')
%irf_legend(h(ip),'T_{\perp}',[0.8 0.1],'color','w')
ylabel(h(ip),{'W_e','[eV]'},'Interpreter','tex');
grid(h(ip), 'off');

% 4.7. electron energy flux
ip = ip + 1;
h(ip) = irf_panel('eEnflux_bg');
c_eval('irf_spectrogram(h(ip), eDist?_bg.deflux.omni.specrec(''energy''), ''log'');', ic);
hold(h(ip), 'on');
c_eval('irf_plot(h(ip), scPot?, ''LineWidth'', 1.3);', ic);
hold(h(ip), 'off');
h(ip).YScale = 'log';
h(ip).YTick = 10.^[1 2 3 4];
caxis(h(ip), [3.5 9]);
irf_legend(h(ip),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
%irf_legend(h(ip),'T_{||}',[0.75 0.1],'color','k')
%irf_legend(h(ip),'T_{\perp}',[0.8 0.1],'color','w')
ylabel(h(ip),{'W_e','[eV]'},'Interpreter','tex');
grid(h(ip), 'off');

% 4.X. global
load('caa/cmap.mat');
c_eval('colormap(h(?), cmap);', 5 : 7);
title(h(1), strcat('MMS', num2str(ic)));
irf_pl_mark(h(1 : npanel1), tint, 'yellow');
irf_pl_mark(h(1 : npanel1), tint.start, 'k', 'Linewidth', 0.6);
irf_pl_mark(h(1 : npanel1), tint.stop, 'k', 'Linewidth', 0.6);
irf_plot_axis_align(h(1 : npanel1))
irf_zoom(h(1 : npanel1),'x',Tint);

% ****************************************
% +++++++++++++++++++++++++++++++ eDist projection +++++++++++++++++++++++++++++++
vlim = 2.0e4;   % x and ylim for electron
projclim = [0 4.8];       % colorbar limit
elevlim = 11.25*1;           % angle over plane to include in slice
strCMap = 'jet';        % colormap


% 5.M1 2D projection; Vx-Vy
ip = ip + 1;
xyz = [dir1; dir2; dir3];
h(ip) = irf_panel('Vx-Vy');
vlabels = {dirstr{1}, dirstr{2}, dirstr{3}};
c_eval('[~, hcb1xy] = mms.plot_projection(h(ip), eDist?.convertto(''s^3/km^6''),''tint'',tint, ''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''scpot'',scPot1,''vlabel'',vlabels, ''clim'', projclim);', ic);
hold(h(ip), 'on');
plot(h(ip), [0, hatB1(1)*15], [0, hatB1(2)*15], 'k', 'LineWidth', 2);
%plot(h(ip), gseVExB1tint(1)/1e2, gseVExB1tint(2)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
%plot(h(ip), gseVe1perptmp(1)/1e2, gseVe1perptmp(2)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%plot(h(ip), [0, M0(3)*10], [0, M0(1)*10], 'b', 'LineWidth', 2);
%plot(h(ip), [0, N0(3)*10], [0, N0(1)*10], 'r', 'LineWidth', 2);
hold(h(ip), 'off');
colormap(h(ip),cmap);
delete(hcb1xy);
strtmp = h(ip).Title.String{1};
%h(ip).Title.String = ['MMS ', strtmp(12:23)];
%delete(h(ip).Title);
h(ip).FontSize = 12;
set(h(ip), 'position', [xpos2 ypos2 xwidth2 ywidth2]);
h(ip).Title.String = '';

% 5.M1 2D projection; Vy_Vz
ip = ip + 1;
xyz = [dir2; dir3; dir1];
h(ip) = irf_panel('Vy-Vz');
vlabels = {dirstr{2}, dirstr{3}, dirstr{1}};
c_eval('[~, hcb1yz] = mms.plot_projection(h(ip), eDist?.convertto(''s^3/km^6''),''tint'',tint, ''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''scpot'',scPot1,''vlabel'',vlabels, ''clim'', projclim);', ic);
hold(h(ip), 'on');
plot(h(ip), [0, hatB1(2)*15], [0, hatB1(3)*15], 'k', 'LineWidth', 2);
%plot(h(ip), gseVExB1tint(2)/1e2, gseVExB1tint(3)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
%plot(h(ip), gseVe1perptmp(2)/1e2, gseVe1perptmp(3)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%plot(h(ip), [0, M0(3)*10], [0, M0(1)*10], 'b', 'LineWidth', 2);
%plot(h(ip), [0, N0(3)*10], [0, N0(1)*10], 'r', 'LineWidth', 2);
hold(h(ip), 'off');
colormap(h(ip),cmap);
delete(hcb1yz);
strtmp = h(ip).Title.String{1};
%h(ip).Title.String = ['MMS ', strtmp(12:23)];
%delete(h(ip).Title);
h(ip).FontSize = 12;
set(h(ip), 'position', [xpos2+xwidth2-0.076 ypos2 xwidth2 ywidth2]);
h(ip).Title.String = '';

% 5.M1 2D projection; Vx_Vz
ip = ip + 1;
xyz = [dir1; dir3; dir2*(-1)];
h(ip) = irf_panel('Vx-Vz');
vlabels = {dirstr{1}, dirstr{3}, dirstr{2}};
c_eval('[~, hcb1xz] = mms.plot_projection(h(ip), eDist?.convertto(''s^3/km^6''),''tint'',tint, ''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''scpot'',scPot1,''vlabel'',vlabels, ''clim'', projclim);', ic);
hold(h(ip), 'on');
plot(h(ip), [0, hatB1(1)*15], [0, hatB1(3)*15], 'color', 'k', 'LineWidth', 2);
%plot(h(ip), gseVExB1tint(1)/1e2, gseVExB1tint(3)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
%plot(h(ip), gseVe1perptmp(1)/1e2, gseVe1perptmp(3)/1e2, 'wo', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%plot(h(ip), [0, L0(3)*10], [0, L0(1)*10], 'k', 'LineWidth', 2);
%plot(h(ip), [0, M0(3)*10], [0, M0(1)*10], 'b', 'LineWidth', 2);
%plot(h(ip), [0, N0(3)*10], [0, N0(1)*10], 'r', 'LineWidth', 2);
hold(h(ip), 'off');
hcb1xz.Position = [0.960 0.7299 0.0100 0.2504];
colormap(h(ip),cmap);
strtmp = h(ip).Title.String{1};
%h(ip).Title.String = ['MMS ', strtmp(12:23)];
%delete(h(ip).Title);
h(ip).FontSize = 12;
set(h(ip), 'position', [xpos2+(xwidth2-0.076)*2 ypos2 xwidth2 ywidth2]);
h(ip).Title.String = '';

% 6. Profiles
ip = ip + 1;
h(ip) = irf_panel('PSD profile');
hplot1 = semilogy(h(ip), Venergy_tint/1e6, eDist_para_tmp, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
hold(h(ip), 'on');
hplot2 = semilogy(h(ip), Venergy_tint/1e6, eDist_perp_tmp, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hplot3 = semilogy(h(ip), Venergy_tint/1e6, eDist_apar_tmp, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
hplot4 = semilogy(h(ip), Venergy_tint/1e6, eDist_bg_para_tmp, 'kv-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
hplot5 = semilogy(h(ip), Venergy_tint/1e6, eDist_bg_perp_tmp, 'bv-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hplot6 = semilogy(h(ip), Venergy_tint/1e6, eDist_bg_apar_tmp, 'rv-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
%hplot7 = semilogy(h(ip), Venergy_tint/1e6, eDist_bgrm_para_tmp, 'k^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
%hplot8 = semilogy(h(ip), Venergy_tint/1e6, eDist_bgrm_perp_tmp, 'b^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
%hplot9 = semilogy(h(ip), Venergy_tint/1e6, eDist_bgrm_apar_tmp, 'r^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

%c_eval('hplot? = semilogy(h(ip), eDist_pad_pa, eDist_pad_tint(?, :)*1e30, ''o-'', ''LineWidth'', 1.5);', 2:26);
hold(h(ip), 'off');
%h(ip).XLim = [0 180];
h(ip).YLim = [5e-2 9e4];
h(ip).XLim = [0 14];
h(ip).XLabel.String = 'V_e [km/s]';
h(ip).YLabel.String = 'f_e [s^3/km^6]';
%h(ip).YTick = '';
set(h(ip), 'position', [xpos3 ypos3 xwidth3+0.05 ywidth3]);
ylabel_step = 0.04;
irf_legend(h(ip), 'f_{||}', [0.71 0.95], 'color', 'k', 'fontsize', 16, 'fontweight', 'bold');
irf_legend(h(ip), 'f_{\perp}', [0.77 0.95], 'color', 'b', 'fontsize', 16, 'fontweight', 'bold');
irf_legend(h(ip), 'f_{-||}', [0.83 0.95], 'color', 'r', 'fontsize', 16, 'fontweight', 'bold');
legend(h(ip), [hplot1 hplot4], {'f_{local}', 'f_{photo}'});

%c_eval('irf_legend(h(ip), energy_str{?}, [1.01 1.00-ylabel_step*?], ''color'', hplot?.Color, ''fontsize'', 10, ''fontWeight'', ''bold'');', 1:23);
%%