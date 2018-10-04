% A routine to compute and plot reduced ion distributions from FPI
%
% The example is fairly slow. Approx 2 min. The code will slow down
% significantly when it reaches the magnetosheath.
%
% Written by A. Johlander


%% Set parameters and get data
% time interval
tint = irf.tint('2015-12-28T03:57:10/2015-12-28T03:59:00');
t = irf.time_array('2015-12-28T03:57:40.3'); % time for 2D distribution

% sc number
ic = 2;

% define velocity grid
vg1D = linspace(-800,800,100);
vg2D = linspace(-1500,1500,200);

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 2e2;

% velocity limit in plot
vlim = 800; % km/s

% get ion distribution (mms.get_data is somewhat slow)
% also get errors
iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);

% ignore flux where count is 1 (also makes function faster)
iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;

% get magnetic field in dmpa
B = mms.get_data('B_dmpa_fgm_brst_l2',tint,ic);

% shock normal vector in GSE 
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
c_eval('nDMPA = mms_dsl2gse(nGSE,defatt?);', ic);
c_eval('t1DMPA = mms_dsl2gse(t1GSE,defatt?);', ic);
c_eval('t2DMPA = mms_dsl2gse(t2GSE,defatt?);', ic);


%% Reduce distribution along normal vector
tic
% reduced distribution along normal vector
f1Dn = iPDist.reduce('1D',nDMPA,'vg',vg1D,'nMC',nMC); 
toc

%% Reduce ion distributions in 2d planes (n-t1 & n-t2)
indT = interp1(iPDist.time.epochUnix,1:iPDist.length,t.epochUnix,'nearest');
% take average for 5 time steps
f2Dnt1 = iPDist(indT-2:indT+2).reduce('2D',nDMPA,t1DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);
f2Dnt2 = iPDist(indT-2:indT+2).reduce('2D',nDMPA,t2DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);
f2Dt1t2 = iPDist(indT-2:indT+2).reduce('2D',t1DMPA,t2DMPA,'base','cart','vg',vg2D,'nMC',nMC*5);


%% Plot reduced distribution as a time series
% make figure
h = irf_plot(2,'newfigure');

% fix figure
% -----------------------------------------------------
set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 9;
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
irf_plot(hca,B)
hold(hca,'on')
irf_plot(hca,B.abs)
hca.YLimMode = 'auto';
ylabel(hca,'$B_{dmpa}$ [nT]','interpreter','latex')
irf_legend(hca,{'B_x';'B_y';'B_z';'|B|'},[1.02,0.9])

% Plot reduced distribution
hca = irf_panel(h,'pdist');
irf_spectrogram(hca,f1Dn.specrec('1D_velocity'),'donotshowcolorbar');
hold(hca,'on')
hcb = colorbar(hca);
ylabel(hcb,'$\log_{10} F_i$ [s m$^{-4}$]','interpreter','latex','fontsize',15)
ylabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
colormap(hca,'jet')
irf_zoom(hca,'y',[min(vg1D),max(vg1D)])


% More figure things
irf_plot_axis_align(h)
irf_zoom(h,'x',tint)
for ii = 1:length(h)
    h(ii).FontSize = 15;
    h(ii).LineWidth = 1.3;
    h(ii).Layer = 'top';
    h(ii).Position(3) = 0.71;
    irf_pl_mark(h(ii),t,'k')
end
h(end).XLabel.String = '';
hcb.LineWidth = 1.3;
hcb.Position(1) = 0.9;



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
xlabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
ylabel(hca,'$V_{t1}$ [km/s]','interpreter','latex')

hca = subplot(1,3,2);
h(2) = hca;
hold(hca,'on')
f2Dnt2.plot_plane(hca,'docolorbar',0)
axis(hca,'equal')
hca.XLim = 1000*[-1,1];
hca.YLim = 1000*[-.5,1.5];
xlabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
ylabel(hca,'$V_{t2}$ [km/s]','interpreter','latex')

hca = subplot(1,3,3);
h(3) = hca;
hold(hca,'on')
f2Dt1t2.plot_plane(hca,'docolorbar',0)
axis(hca,'equal')
hca.XLim = 1000*[-1,1];
hca.YLim = 1000*[-.5,1.5];
xlabel(hca,'$V_{t1}$ [km/s]','interpreter','latex')
ylabel(hca,'$V_{t2}$ [km/s]','interpreter','latex')

% colorbar
htempPos = h(1).Position;
hcb = colorbar(h(1),'north');
h(1).Position = htempPos;
hcb.Position = [h(1).Position(1),0.83,sum(h(end).Position([1,3]))-h(1).Position(1),0.04];
ylabel(hcb,'$\log_{10} F_i$ [s$^2$ m$^{-5}$]','interpreter','latex','fontsize',15)
hcb.LineWidth = 1.3;

for ii = 1:length(h)
    % plot crosseye at origin in spacecrat frame
    plot(h(ii),[h(ii).XLim,0,0,0],[0,0,0,h(ii).YLim],'k-','Linewidth',1.2)
    h(ii).LineWidth = 1.3;
    h(ii).Layer = 'top';
    % match clim
    h(ii).CLim = h(1).CLim;
    h(ii).Position(2) = 0.06;
end

