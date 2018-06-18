% A routine to compute and plot reduced electron distributions from FPI
% 
% Compare with: Wilder, F. D., et al. (2016), GRL, 43, 5909?5917,
% doi:10.1002/2016GL069473.
%
% The Example is fairly slow. Approx 4 min.
%
% Written by A. Johlander


%% Set parameters and get data
% time interval
tint = irf.tint('2015-09-19T10:08:14/2015-09-19T10:08:18');
% times to make lines
t1 = irf.time_array('2015-09-19T10:08:16.275');
t2 = irf.time_array('2015-09-19T10:08:16.335');

% sc number
ic = 4;

% color/y-limit
clim = 10.^[-4.5,-0.5]; % s m^-4

% define velocity grid
vg = linspace(-45e3,45e3,100); % km/s  

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 1e3;

% velocity limit in plot
vlim = 50e3; % km/s

% get distribution function
c_eval('ePDist = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint));',ic)
ePDist = ePDist.tlim(tint);

% get magnetic field in DMPA (since ePDist is in DMPA)
c_eval('B = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);',ic)

% Get the spacecraft potential (in this example, it is not so important)
c_eval('scPot = mms.get_data(''V_edp_brst_l2'',tint,?);',ic)

% remove flux from bottom two energy levels to make it more like in the
% paper
ePDist.data(:,1:2,:,:) = 0;


%% Reduce distribution
tic
% reduced distribution along B
f1D = ePDist.reduce('1D',B,'vg',vg,'nMC',nMC,'scpot',scPot); 
toc

%% Plot reduced distribution as a time series
% make figure
h = irf_plot(2,'newfigure');

% plot magnetic field
hca = irf_panel(h,'Bxyz');
irf_plot(hca,B)
ylabel(hca,'Bdmpa [nT]')
irf_legend(hca,{'Bx';'By';'Bz'},[1.02,0.9])

% Plot reduced distribution
hca = irf_panel(h,'pdist');
[~,hcb] = irf_spectrogram(hca,f1D.specrec('1D_velocity'));
hcb.Label.String = 'log_{10}F_e [s m^{-4}]';
ylabel(hca,'V_{||} [km/s]')
colormap('jet')
irf_zoom(hca,'y',[min(vg),max(vg)])
irf_plot_axis_align(h)
irf_zoom(h,'x',tint)


%% plot distribution as lines for the two selected lines

% get indicies for the times
it1 = interp1(ePDist.time.epochUnix,1:length(ePDist),t1.epochUnix,'nearest');
it2 = interp1(ePDist.time.epochUnix,1:length(ePDist),t2.epochUnix,'nearest');

% matlab colors
col = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980];

% initiate figure
hca = irf_plot(1,'newfigure');
hold(hca,'on')
% plot the lines
plot(hca,f1D(it1).depend{1},f1D(it1).data,'Color',col(2,:),'linewidth',2)
plot(hca,f1D(it2).depend{1},f1D(it2).data,'Color',col(1,:),'linewidth',2)
hca.YScale = 'log';
hca.YLim = clim;
% legends show time centers
hca.ColorOrder = flipud(col(1:2,:)); % set colr order for legends
irf_legend(hca,{ePDist(it1).time.toUtc;ePDist(it2).time.toUtc;},[0.98,0.98])
hca.XLim = [min(vg),max(vg)];
%labels
xlabel(hca,'V_{||} [km/s]')
ylabel(hca,'F [s m^-^4]')

