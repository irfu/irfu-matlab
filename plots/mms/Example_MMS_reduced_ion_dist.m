% A routine to compute and plot reduced ion distributions from FPI
%
% The example is fairly slow. Approx 2 min. The code will slow down
% significantly when it reaches the magnetosheath.
%
% Written by A. Johlander


%% Set parameters and get data
% time interval
tint = irf.tint('2015-12-28T03:57:1/2015-12-28T03:59:00');


% sc number
ic = 4;

% define velocity grid
vg = linspace(-800,800,100);

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 2e2;

% velocity limit in plot
vlim = 800; % km/s

% make two PDist objects to get longer time
% also get errors
c_eval('[iPDistA,iPDistAerr] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint(1)));',ic)
c_eval('[iPDistB,iPDistBerr] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint(2)));',ic)

% ignore flux where count is 1 (also makes function faster)
iPDistA.data(iPDistA.data<1.1*iPDistAerr.data) = 0;
iPDistB.data(iPDistB.data<1.1*iPDistBerr.data) = 0;

% get magnetic field in dmpa
c_eval('B = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);',ic)

% shock normal vector in GSE 
nvec = [0.9580   -0.2708   -0.0938]; nvec = nvec/norm(nvec);

% make nvec a TSeries object resampled to B
nGSE = irf.ts_vec_xyz(B.time,repmat(nvec,length(B),1));

% rotate normal vector to DMPA
c_eval('defatt?=mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
c_eval('defatt?.zdec=mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
c_eval('nDMPA = mms_dsl2gse(nGSE,defatt?);', ic);


%% Reduce distribution
tic
% reduced distribution along B
f1DA = iPDistA.reduce('1D',nDMPA,'vg',vg,'nMC',nMC); 
f1DB = iPDistB.reduce('1D',nDMPA,'vg',vg,'nMC',nMC); 
toc

%% Plot reduced distribution as a time series
% make figure
h = irf_plot(2,'newfigure');

% fix figure
% -----------------------------------------------------
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 7;
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
irf_spectrogram(hca,f1DA.specrec,'donotshowcolorbar');
hold(hca,'on')
irf_spectrogram(hca,f1DB.specrec,'donotshowcolorbar');
hcb = colorbar(hca);
ylabel(hcb,'$\log_{10} F_i$ [s m$^{-4}$]','interpreter','latex')
ylabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
colormap('jet')
irf_zoom(hca,'y',[min(vg),max(vg)])

% More figure things
irf_plot_axis_align(h)
irf_zoom(h,'x',[iPDistA.time(1),iPDistB.time(end)])
for ii = 1:length(h)
    h(ii).FontSize = 15;
    h(ii).LineWidth = 1.3;
    h(ii).Layer = 'top';
    h(ii).Position(3) = 0.71;
end
h(end).XLabel.String = '';
hcb.LineWidth = 1.3;
hcb.Position(1) = 0.9;