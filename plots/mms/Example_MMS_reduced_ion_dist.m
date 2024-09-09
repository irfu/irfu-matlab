% A routine to compute and plot reduced ion distributions from FPI
%
% Since this example was written, the event it uses has been published in
% https://doi.org/10.3847/1538-4357/abcb88
%
% See Example_MMS_IPshock for how to convert these type of data to the
% normal incidence frame
%
% Written by A. Johlander


%% Set parameters and get data
% time interval
tint = irf.tint('2015-12-28T03:57:10/2015-12-28T03:59:00');
t = irf.time_array('2015-12-28T03:57:40.3'); % time for 2D distribution

% sc number
ic = 2;

% define velocity grid
vnlim = [-800,800];
vt1lim = vnlim;
vt2lim = vnlim+300;

vg1Dn = linspace(vnlim(1),vnlim(2),100);
vg1Dt1 = linspace(vt1lim(1),vt1lim(2),100);
vg1Dt2 = linspace(vt2lim(1),vt2lim(2),100);
% larger v-space for 2D distribution
vg2D = linspace(-1500,1500,200);

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 2e2;

%% Get data
% get ion distribution (mms.get_data is somewhat slow)
% also get errors
iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);

% ignore psd where count is 1 (also makes function faster)
iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;

% get magnetic field in dmpa
B = mms.get_data('B_dmpa_fgm_brst_l2',tint,ic);

% shock normal vector in GSE (get this from irf_shock_normal or irf_shock_gui)
nvec = [0.9580   -0.2708   -0.0938]; nvec = nvec/norm(nvec);
% Upstream magnetic field in GSE
Bu = [-1.0948   -2.6270    1.6478];
% t2 vector in GSE (same vectors as in [Johlander et al. 2016, PRL])
t2vec = cross(nvec,Bu)/norm(cross(nvec,Bu));
% t1 vector in GSE
t1vec = cross(t2vec,nvec);


% make vectors to TSeries objects resampled to B
nGSE = irf.ts_vec_xyz(B.time,repmat(nvec,length(B),1));
t1GSE = irf.ts_vec_xyz(B.time,repmat(t1vec,length(B),1));
t2GSE = irf.ts_vec_xyz(B.time,repmat(t2vec,length(B),1));

% rotate vectors to DMPA
c_eval('defatt?=mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
c_eval('defatt?.zdec=mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
c_eval('nDMPA = mms_dsl2gse(nGSE,defatt?,-1);', ic);
c_eval('t1DMPA = mms_dsl2gse(t1GSE,defatt?,-1);', ic);
c_eval('t2DMPA = mms_dsl2gse(t2GSE,defatt?,-1);', ic);


%% Reduce distribution along the three vectors
tic
% reduced distributions
f1Dn = iPDist.reduce('1D',nDMPA,'vg',vg1Dn,'nMC',nMC);
f1Dt1 = iPDist.reduce('1D',t1DMPA,'vg',vg1Dt1,'nMC',nMC);
f1Dt2 = iPDist.reduce('1D',t2DMPA,'vg',vg1Dt2,'nMC',nMC);
toc

%% Reduce ion distributions in 2d planes (n-t1 & n-t2)
indT = interp1(iPDist.time.epochUnix,1:iPDist.length,t.epochUnix,'nearest');
% take average for 5 time steps
f2Dnt1 = iPDist(indT-2:indT+2).reduce('2D',nDMPA,t1DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);
f2Dnt2 = iPDist(indT-2:indT+2).reduce('2D',nDMPA,t2DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);
f2Dt1t2 = iPDist(indT-2:indT+2).reduce('2D',t1DMPA,t2DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);


%% set colormap
% nice colormap but opens a figure for some reason
% cmap = irf_colormap('waterfall');
% just jet
cmap = jet;


%% Plot reduced distribution as a time series
% make figure
h = irf_plot(4,'newfigure');

% fix figure
% -----------------------------------------------------
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 14;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
% additional good options
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
% -----------------------------------------------------

% plot magnetic field
hca = irf_panel(h,'Bxyz');
irf_plot(hca,B,'linewidth',1.5)
hold(hca,'on')
irf_plot(hca,B.abs,'linewidth',1.5)
hlinetemp = irf_plot(hca,[tint.epochUnix,[0;0]],'--','color',[1,1,1]*.5,'linewidth',1.5);
uistack(hlinetemp,'bottom') % put it at the bottom
hca.YLimMode = 'auto';
ylabel(hca,'$\mathbf{B}_{\mathrm{dmpa}}$ [nT]','interpreter','latex')
irf_legend(hca,{'$B_x$';'$B_y$';'$B_z$';'$|\mathbf{B}|$'},[1.02,0.95],...
  'interpreter','latex','fontsize',18)

% Plot reduced distribution along n
hca = irf_panel(h,'f1dn');
irf_spectrogram(hca,f1Dn.specrec('1D_velocity'),'donotshowcolorbar');
hold(hca,'on')
ylabel(hca,'$v_{n}$ [km\,s$^{-1}$]','interpreter','latex')
colormap(hca,'jet')
irf_zoom(hca,'y',vnlim)

% Plot reduced distribution along t1
hca = irf_panel(h,'f1dt1');
irf_spectrogram(hca,f1Dt1.specrec('1D_velocity'),'donotshowcolorbar');
hold(hca,'on')
ylabel(hca,'$v_{t1}$ [km\,s$^{-1}$]','interpreter','latex')
colormap(hca,'jet')
irf_zoom(hca,'y',vt1lim)

% Plot reduced distribution along t2
hca = irf_panel(h,'f1dt2');
irf_spectrogram(hca,f1Dt2.specrec('1D_velocity'),'donotshowcolorbar');
hold(hca,'on')
ylabel(hca,'$v_{t2}$ [km\,s$^{-1}$]','interpreter','latex')
colormap(hca,'jet')
irf_zoom(hca,'y',vt2lim)


% More figure things

irf_pl_number_subplots(h,[.03,.93],'num','\bf{(?)}','fontsize',18,...
  'interpreter','latex','backgroundcolor','w','edgecolor','k','linewidth',1.8)

irf_plot_axis_align(h)
irf_zoom(h,'x',tint)
for ii = 1:length(h)
  h(ii).FontSize = 18;
  h(ii).LineWidth = 1.8;
  h(ii).Layer = 'top';
  h(ii).Position(1) = 0.135;
  h(ii).Position(3) = 0.71;
  irf_pl_mark(h(ii),t,'k','linewidth',1.5)
  grid(h(ii),'off')
  h(ii).YLabel.Units = 'normalized';
  h(ii).YLabel.Position(1) = -.125;
  colormap(h(ii),cmap)
end

%
% ---- colorbar ----
% first fix the color limit for all panels
h(3).CLim = h(2).CLim;
h(4).CLim = h(2).CLim;

hcb = colorbar(h(3));
ylabel(hcb,'$\log_{10} F_i$ [s\,m$^{-4}$]','interpreter','latex','fontsize',18)
hcb.LineWidth = 1.8;
hcb.Position([2,4]) = [h(4).Position(2),3*h(4).Position(4)];
hcb.Position(1) = sum(h(4).Position([1,3]))+0.015;

h(end).XTickLabelRotation = 0;


%% Plot 2D distributions
% make figure
fig = figure;
fig.Color = 'w';
h = gobjects(1,2);

% fix figure
% -----------------------------------------------------
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
% additional good options
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
% -----------------------------------------------------


hca = subplot(1,3,1);
h(1) = hca;
hold(hca,'on')
f2Dnt1.plot_plane(hca,'docolorbar',0)
axis(hca,'equal')
hca.XLim = 1000*[-1,1];
hca.YLim = 1000*[-1,1];
xlabel(hca,'$v_{n}$ [km\,s$^{-1}$]','interpreter','latex')
ylabel(hca,'$v_{t1}$ [km\,s$^{-1}$]','interpreter','latex')

hca = subplot(1,3,2);
h(2) = hca;
hold(hca,'on')
f2Dnt2.plot_plane(hca,'docolorbar',0)
axis(hca,'equal')
hca.XLim = 1000*[-1,1];
hca.YLim = 1000*[-.5,1.5];
xlabel(hca,'$v_{n}$ [km\,s$^{-1}$]','interpreter','latex')
ylabel(hca,'$v_{t2}$ [km\,s$^{-1}$]','interpreter','latex')

hca = subplot(1,3,3);
h(3) = hca;
hold(hca,'on')
f2Dt1t2.plot_plane(hca,'docolorbar',0)
axis(hca,'equal')
hca.XLim = 1000*[-1,1];
hca.YLim = 1000*[-.5,1.5];
xlabel(hca,'$v_{t1}$ [km\,s$^{-1}$]','interpreter','latex')
ylabel(hca,'$v_{t2}$ [km\,s$^{-1}$]','interpreter','latex')

% ---- colorbar ----
% first fix the color limit for all panels
h(2).CLim = h(1).CLim;
h(3).CLim = h(1).CLim;
htempPos = h(1).Position;
hcb = colorbar(h(1),'north');
h(1).Position = htempPos;
hcb.Position = [h(1).Position(1),0.80,sum(h(end).Position([1,3]))-h(1).Position(1),0.04];
ylabel(hcb,'$\log_{10} F_i$ [s$^2$\,m$^{-5}$]','interpreter','latex','fontsize',18)
hcb.LineWidth = 1.8;

for ii = 1:length(h)
  % plot crosseye at origin in spacecrat frame
  plot(h(ii),[h(ii).XLim,0,0,0],[0,0,0,h(ii).YLim],'k-','Linewidth',1.3)
  h(ii).LineWidth = 1.8;
  h(ii).FontSize = 18;
  h(ii).Layer = 'top';
  colormap(h(ii),cmap)
  % match clim
  h(ii).CLim = h(1).CLim;
  h(ii).Position(2) = 0.06;
end

h(1).Position(1) = .1;
h(3).Position(1) = .72;

