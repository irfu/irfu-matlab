% A routine to compute and plot reduced electron distributions from FPI
%
% The function computes 2D reduced distributions from all four MMS. The
% plotted distribution function is averaged both over time and spacecraft.
%
% The Example is fairly slow. 
%
% Written by A. Johlander


%% Set parameters and get data
% time interval (0.5 s)
tint = irf.tint('2017-07-26T07:01:46.25/2017-07-26T07:01:46.75');

% define velocity grid
vg = linspace(-70e3,70e3,100); % km/s  

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 2e3;

% It is possible to limit the out-of-plane speed, in a way making a cut
% (This seems to be undocumented)
vzint = 20e3*[-1,1]; % km/s

% get distribution function for all four MMS
c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint));')
c_eval('ePDist? = ePDist?.tlim(tint);')

% get magnetic ield in DMPA (since ePDist is in DMPA)
c_eval('B? = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);')

% get average magnetic field
B = (B1+B2.resample(B1)+B3.resample(B1)+B4.resample(B1))/4;

% use a fixed coordinate system
% unit vector pointing along B
bvec = mean(B.tlim(tint).data)/norm(mean(B.tlim(tint).data));
% perp unit vectors, (par,perp1,perp2) is right handed
perp1vec = cross([0,0,1],bvec); perp1vec = perp1vec/norm(perp1vec);
perp2vec = cross(bvec,perp1vec);


% remove flux from bottom 11 energy levels to get rid of internal
% photoelectrons
c_eval('ePDist?.data(:,1:11,:,:) = 0;')


%% Reduce distributions 
tic
c_eval('f2Dparperp1? = ePDist?.reduce(''2D'',bvec,perp1vec,''base'',''cart'',''vg'',vg,''nMC'',nMC,''vint'',vzint);')
c_eval('f2Dparperp2? = ePDist?.reduce(''2D'',bvec,perp2vec,''base'',''cart'',''vg'',vg,''nMC'',nMC,''vint'',vzint);')
c_eval('f2Dperp1perp2? = ePDist?.reduce(''2D'',perp1vec,perp2vec,''base'',''cart'',''vg'',vg,''nMC'',nMC,''vint'',vzint);')
toc


%% get spacecraft-averaged reduced distriubtions (pretty much a hack right now)
% code should be improved

% PDist object with reduced distribution averaged over time and spacecraft
tempDistData = (mean(f2Dparperp11.data)+mean(f2Dparperp12.data)+mean(f2Dparperp13.data)+mean(f2Dparperp14.data))/4;
f2Dparperp1 = PDist(tint(1),tempDistData,'plane (reduced)',f2Dparperp11.depend{1}(1,:),f2Dparperp11.depend{2}(1,:));

tempDistData = (mean(f2Dparperp21.data)+mean(f2Dparperp22.data)+mean(f2Dparperp23.data)+mean(f2Dparperp24.data))/4;
f2Dparperp2 = PDist(tint(1),tempDistData,'plane (reduced)',f2Dparperp21.depend{1}(1,:),f2Dparperp21.depend{2}(1,:));

tempDistData = (mean(f2Dperp1perp21.data)+mean(f2Dperp1perp22.data)+mean(f2Dperp1perp23.data)+mean(f2Dperp1perp24.data))/4;
f2Dperp1perp2 = PDist(tint(1),tempDistData,'plane (reduced)',f2Dperp1perp21.depend{1}(1,:),f2Dperp1perp21.depend{2}(1,:));

% set species and ancillary
f2Dparperp1.species = 'electrons';
f2Dparperp1.ancillary = f2Dparperp11.ancillary;

f2Dparperp2.species = 'electrons';
f2Dparperp2.ancillary = f2Dparperp21.ancillary;

f2Dperp1perp2.species = 'electrons';
f2Dperp1perp2.ancillary = f2Dperp1perp21.ancillary;


%% Plot 2D distributions of single, time-averaged, and time-spacecraft-averaged
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
f2Dparperp11(1).plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\parallel}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,1}$ [$10^3$km/s]','interpreter','latex')
title(hca,'One distribution');

hca = subplot(1,3,2);
h(2) = hca;
f2Dparperp11.plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\parallel}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,1}$ [$10^3$km/s]','interpreter','latex')
title(hca,'Time-averaged, MMS1');

hca = subplot(1,3,3);
h(3) = hca;
f2Dparperp1.plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\parallel}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,1}$ [$10^3$km/s]','interpreter','latex')
title(hca,'Time- and spacecraft-averaged');

% colorbar
htempPos = h(1).Position;
hcb = colorbar(h(1),'north');
h(1).Position = htempPos;
hcb.Position = [h(1).Position(1),0.83,sum(h(end).Position([1,3]))-h(1).Position(1),0.04];
ylabel(hcb,'$\log_{10} F_e$ [s$^2$ m$^{-5}$]','interpreter','latex','fontsize',15)
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

% new colormap!
irf_colormap('waterfall')

%% Plot time and spacraft averaged 2D distributions
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
f2Dparperp1.plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\parallel}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,1}$ [$10^3$km/s]','interpreter','latex')

hca = subplot(1,3,2);
h(2) = hca;
f2Dparperp2.plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\parallel}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,2}$ [$10^3$km/s]','interpreter','latex')

hca = subplot(1,3,3);
h(3) = hca;
f2Dperp1perp2.plot_plane(hca,'colorbar',0,'contour',6)
hold(hca,'on')
axis(hca,'equal')
xlabel(hca,'$V_{\perp,1}$ [$10^3$km/s]','interpreter','latex')
ylabel(hca,'$V_{\perp,2}$ [$10^3$km/s]','interpreter','latex')

% colorbar
htempPos = h(1).Position;
hcb = colorbar(h(1),'north');
h(1).Position = htempPos;
hcb.Position = [h(1).Position(1),0.76,sum(h(end).Position([1,3]))-h(1).Position(1),0.04];
ylabel(hcb,'$\log_{10} F_e$ [s$^2$ m$^{-5}$]','interpreter','latex','fontsize',15)
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

ht = title(h(2),'Four-spacecraft averages');
ht.Units = 'normalized';
ht.Position(2) = 1.35;

% new colormap!
irf_colormap('waterfall')
