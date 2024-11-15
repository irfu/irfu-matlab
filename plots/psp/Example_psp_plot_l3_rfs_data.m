% Plots Level 3 data from PSP's Radio Frequency Spectrometer (RFS)
% Both the High frequency receiver (HFR) and Low frequency receiver (LFR) bands

%% Choosing the day
% clear all
% % Set 'global' if you want to load data without declaring the dateStart and dateEnd args
% global dateStart dateEnd

% Tint = irf.tint('2020-06-07T00:00:00.000Z/2020-06-07T23:59:59.999Z') %#ok<NOPTS>
Tint = irf.tint('2024-09-09T00:00:00.000Z/2024-09-09T23:59:59.999Z');
time_str = toUtc(Tint,1);
dateStart = [time_str(1,1:4) ' ' time_str(1,6:7) ' ' time_str(1,9:10)]; % String vector
dateEnd = [time_str(2,1:4) ' ' time_str(2,6:7) ' ' time_str(2,9:10)];


psp_load([],'l3_rfs_hfr',dateStart,dateEnd) % uses data store
% psp_load('./','l3_rfs_hfr',dateStart,dateEnd) % If the data is in the same folder

rfs_hfr_v1v2_spec = struct('t', l3_rfs_hfr_v1v2_freq.time.epochUnix);
rfs_hfr_v1v2_spec.p = l3_rfs_hfr_v1v2.data;
rfs_hfr_v1v2_spec.f = single(l3_rfs_hfr_v1v2_freq.data);
rfs_hfr_v1v2_spec.p_label = [ '[ log_{10}' '(V^2/Hz)' ']' ];
rfs_hfr_v1v2_spec.f_label = 'f [Hz]';

rfs_hfr_v1v2_flux = struct('t', l3_rfs_hfr_v1v2_freq.time.epochUnix);
rfs_hfr_v1v2_flux.p = l3_rfs_hfr_v1v2.data;
rfs_hfr_v1v2_flux.f = single(l3_rfs_hfr_v1v2_freq.data);
rfs_hfr_v1v2_flux.p_label = [ '[ log_{10}' '(W/m^{-2}Hz^{-1})' ']' ];
rfs_hfr_v1v2_flux.f_label = 'f [Hz]';

rfs_hfr_v3v4_spec = struct('t', l3_rfs_hfr_v3v4_freq.time.epochUnix);
rfs_hfr_v3v4_spec.p = l3_rfs_hfr_v3v4.data;
rfs_hfr_v3v4_spec.f = single(l3_rfs_hfr_v3v4_freq.data);
rfs_hfr_v3v4_spec.p_label = [ '[ log_{10}' '(V^2/Hz)' ']' ];
rfs_hfr_v3v4_spec.f_label = 'f [Hz]';

rfs_hfr_v3v4_flux = struct('t', l3_rfs_hfr_v3v4_freq.time.epochUnix);
rfs_hfr_v3v4_flux.p = l3_rfs_hfr_v3v4.data;
rfs_hfr_v3v4_flux.f = single(l3_rfs_hfr_v3v4_freq.data);
rfs_hfr_v3v4_flux.p_label = [ '[ log_{10}' '(W/m^{-2}Hz^{-1})' ']' ];
rfs_hfr_v3v4_flux.f_label = 'f [Hz]';

psp_load([],'l3_rfs_lfr',dateStart,dateEnd)
% psp_load('./','l3_rfs_lfr',dateStart,dateEnd)

rfs_lfr_v1v2_spec = struct('t', l3_rfs_lfr_v1v2_freq.time.epochUnix);
rfs_lfr_v1v2_spec.p = l3_rfs_lfr_v1v2.data;
rfs_lfr_v1v2_spec.f = single(l3_rfs_lfr_v1v2_freq.data);
rfs_lfr_v1v2_spec.p_label = [ '[ log_{10}' '(V^2/Hz)' ']' ];
rfs_lfr_v1v2_spec.f_label = 'f [Hz]';

rfs_lfr_v1v2_flux = struct('t', l3_rfs_lfr_v1v2_freq.time.epochUnix);
rfs_lfr_v1v2_flux.p = l3_rfs_lfr_v1v2.data;
rfs_lfr_v1v2_flux.f = single(l3_rfs_lfr_v1v2_freq.data);
rfs_lfr_v1v2_flux.p_label = [ '[ log_{10}' '(W/m^{-2}Hz^{-1})' ']' ];
rfs_lfr_v1v2_flux.f_label = 'f [Hz]';

rfs_lfr_v3v4_spec = struct('t', l3_rfs_lfr_v3v4_freq.time.epochUnix);
rfs_lfr_v3v4_spec.p = l3_rfs_lfr_v3v4.data;
rfs_lfr_v3v4_spec.f = single(l3_rfs_lfr_v3v4_freq.data);
rfs_lfr_v3v4_spec.p_label = [ '[ log_{10}' '(V^2/Hz)' ']' ];
rfs_lfr_v3v4_spec.f_label = 'f [Hz]';

rfs_lfr_v3v4_flux = struct('t', l3_rfs_lfr_v3v4_freq.time.epochUnix);
rfs_lfr_v3v4_flux.p = l3_rfs_lfr_v3v4.data;
rfs_lfr_v3v4_flux.f = single(l3_rfs_lfr_v3v4_freq.data);
rfs_lfr_v3v4_flux.p_label = [ '[ log_{10}' '(W/m^{-2}Hz^{-1})' ']' ];
rfs_lfr_v3v4_flux.f_label = 'f [Hz]';

%% Plot

Tintt = Tint;
% Tintt = irf.tint('2024-09-09T08:00:00.000Z/2024-09-09T14:00:00.000Z'); zoom-in in a shorter interval
% irf_zoom(h,'x',Tintt);

txt = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z']';
f = figure;
f.Units = 'pixels';
f.Color = [1 1 1];
f.WindowStyle = 'normal';
f.Resize = 'on';
f.ToolBar = 'none';
%%

% Two monitors position
% To adjust, resize to your prefered size and aspect ratio, do f.get on Matlab comand line,
% copy the new vectors for the parameters below and replace them
f.Position = [425 417 1815 896];
f.InnerPosition = [425 417 1815 896];
f.OuterPosition = [425 417 1815 923];

% Color-blind friendly palette
colorblind = ['#000000' % black
             '#5987b5'  % blue
             '#9e3f4f'  % moderate red
             '#ddaa33'  % orange
             '#D6BA1D'  % sunflower yellow
             '#AA3377'  % dark pink
             '#0077BB'  % strong blue
             '#2e913f'  % green
             '#EE7733'];
colororder(colorblind)

fs = 20;
fss = 18;

h = irf_plot(2);

xstart = 0.0475;
xwidth = 0.925;
ystart = 0.95;
ywidth = 0.425;

set(h(1),'position',[xstart ystart-ywidth xwidth ywidth]);
set(h(2),'position',[xstart ystart-2*ywidth xwidth ywidth]);

hfr = 1;
[hca, hcb] = irf_spectrogram(h(hfr),rfs_hfr_v1v2_spec); % ch0 V1V2
% [hca, hcb] = irf_spectrogram(h(hfr),rfs_hfr_v3v4_spec); % ch1 V3V4
% [hca, hcb] = irf_spectrogram(h(hfr),rfs_hfr_v1v2_flux); % Power flux (W/m^{-2}Hz^{-1})
% [hca, hcb] = irf_spectrogram(h(hfr),rfs_hfr_v3v4_flux);
irf_zoom(h(hfr),'y',[1.275e6 1.92e7])
h(hfr).YLimMode = 'auto';
hca.YScale = 'log';
hca.Layer = 'top';
hca.TickDir = 'out';
hca.TickLength = [0.005 0.01];
hca.XAxis.MinorTick = 'on';
hca.YAxis.MinorTick = 'on';
hca.FontName = 'Times New Roman';
hca.FontSize = fs;
hcb.FontSize = fss;
hcb.Label.FontSize = fss;
hcb.Label.FontName = 'Times New Roman';
hcb.TickDirection = 'in';
% caxis(h(hfr),'auto')
clim(h(hfr),[-16.0, -13.5]);
% set(hcb,'Ticks', [-18 -17 -16 -15 -14 -13 -12])
ylabel(h(hfr),'f [Hz]','Interpreter','tex','FontSize',fs,'FontName','Times New Roman')
irf_legend(h(hfr),['(' txt(hfr) ')'],[0.99 0.97],'color','#ffffff','fontsize',fs,'fontweight','bold','FontName','Times New Roman')

lfr = 2;
[hcc, hcd] = irf_spectrogram(h(lfr),rfs_lfr_v1v2_spec); % ch0 V1V2
% [hcc, hcd] = irf_spectrogram(h(lfr),rfs_lfr_v3v4_spec); % ch1 V3V4
% [hcc, hcd] = irf_spectrogram(h(lfr),rfs_lfr_v2v2_flux); % Power flux (W/m^{-2}Hz^{-1})
% [hcc, hcd] = irf_spectrogram(h(lfr),rfs_lfr_v3v4_flux);
irf_zoom(h(lfr),'y',[4.9e4 1.274e6]) 
h(lfr).YAxis.ExponentMode = 'manual';
h(lfr).YAxis.Exponent = 0;
hcc.YScale = 'log';
hcc.YTick = [1e5 1e6];
hcc.Layer = 'top';
hcc.TickDir = 'out';
hcc.TickLength = [0.005 0.01];
hcc.XAxis.MinorTick = 'on';
hcc.YAxis.MinorTick = 'on';
hcc.FontName = 'Times New Roman';
hcc.FontSize = fs;
hcd.FontSize = fss;
hcd.Label.FontSize = fss;
hcd.Label.FontName = 'Times New Roman';
hcd.TickDirection = 'in';
% clim(h(lfr),'auto')
clim(h(lfr),[-16, -13.5]);
% set(hcd,'Ticks', [-18 -17 -16 -15 -14 -13 -12 -11])
ylabel(h(lfr),{'f [Hz]'},'Interpreter','tex','FontSize',fs,'FontName','Times New Roman')
irf_legend(h(lfr),['(' txt(lfr) ')'],[0.99 0.97],'color','#ffffff','fontsize',fs,'fontweight','bold','FontName','Times New Roman')

grid('off')

clmap = irf_colormap('batlow');
colormap(clmap);
% % The colormaps below are color-blind friendly and require
% % the "MatPlotLib Perceptually Uniform Colormaps" add-on
% colormap('inferno');
% colormap('viridis')
% colormap('cividis')
% colormap('magma')
% colormap('plasma')

set(hcb,'position',get(hcb,'position')+[-0.015 0.0 0.0 0.0]); 
set(hcd,'position',get(hcd,'position')+[-0.015 0.0 0.0 0.0]); 

title(h(1),'PSP FIELDS Radio Frequency Spectrometer (RFS), HFR and LFR V1V2 data.','FontSize',fss,'fontweight','normal')
% title(h(1),'PSP FIELDS Radio Frequency Spectrometer (RFS), HFR and LFR V3V4 data.','FontSize',fss,'fontweight','normal')

irf_zoom(h,'y')
xtickangle(h,0)
irf_zoom(h,'x',Tintt)
set(h,'TickDir','in','TickLength',[0.0075 0.02],'XMinorTick','on','YMinorTick','on','FontSize',fs,'FontName','Times New Roman');
irf_plot_axis_align
irf_plot_ylabels_align

time_str = toUtc(Tintt,1);
fig_name = [time_str(1,1:10) '_rfs_v1v2_spec']; % full day
figname = ['/home/you/working_directory/figure_directory/' fig_name '.png']; % adapt this to your directory
% ax = f; exportgraphics(ax,figname,'Resolution',600); % save figure with minimal canvas padding
