% Calculate and plot electron pitch angle distributions from L1b particle 
% brst data. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number

mms.db_init('local_file_db','/data/mms');

tint = irf.tint('2015-09-08T10:29:20.000000Z/2015-09-08T10:29:40.000000Z');
%%

% Load Magnetic Field
c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);

% Load Ion Moments 
c_eval('ni=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',tint);',ic);
dtmoments = ni.time(2)-ni.time(1);
fsmoments = 1/dtmoments;
ni = irf_filt(ni,0,fsmoments/4,fsmoments,3);
c_eval('Vx=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkX'',tint);',ic);
c_eval('Vy=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkY'',tint);',ic);
c_eval('Vz=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkZ'',tint);',ic);
V = TSeries(Vx.time,[Vx.data Vy.data Vz.data],'to',1);
V = irf_filt(V,0,fsmoments/4,fsmoments,3);

% Load Electron Moments
dsStr = sprintf('mms%d_fpi_brst_l1b_des-moms',ic);
ne=mms.db_get_ts(dsStr,sprintf('mms%d_des_numberDensity',ic),tint);
dtmoments = ne.time(2)-ne.time(1);
fsmoments = 1/dtmoments;
ne = irf_filt(ne,0,fsmoments/4,fsmoments,3);

% Load Electron distribution
tic;
c_eval('dist = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ic);
toc

tic;
[paddist,theta] = mms.get_pitchangledist(dist,Bxyz);
toc;


%Convert to pad deflux
paddist.data = paddist.data*1e30;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;
energyspec = ones(length(paddist.time),1)*energy;

for ii = 1:length(theta);
    paddist.data(:,:,ii) = paddist.data(:,:,ii).*energyspec.^2;
end
paddist.data = paddist.data/1e6/(5.486e-4)^2/0.53707;



DEFluxomni = zeros(length(dist.time),length(energy));
for ii = 1:length(dist.time);
    disttemp = squeeze(dist.data(ii,:,:,:));
    DEFluxomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
end
DEFluxomni = DEFluxomni.*energyspec.^2;
DEFluxomni = DEFluxomni/1e6/(5.486e-4)^2/0.53707*1e30;

specomni=struct('t',dist.time.epochUnix);
specomni.p = DEFluxomni;
specomni.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specomni.f_label={''};
specomni.f = energy;

Ethres = 2000;
Elarge = find(energy > Ethres);
Esmall = find(energy < Ethres);

padlargeE = squeeze(mean(paddist.data(:,Elarge,:),2));
padsmallE = squeeze(mean(paddist.data(:,Esmall,:),2));

specpadl=struct('t',paddist.time.epochUnix);
specpadl.p = padlargeE;
specpadl.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
specpadl.f_label={''};
specpadl.f = theta;

specpads = specpadl;
specpads.p = padsmallE;

h=irf_plot(6,'newfigure');

Bxyz = Bxyz.tlim(tint);
h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Vxyz');
irf_plot(h(2),V);
ylabel(h(2),'V_{DSL} (km s^{-1})','Interpreter','tex');
irf_legend(h(2),{'V_{x}','V_{y}','V_{z}'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('ne');
irf_plot(h(3),ne);
hold(h(3),'on');
irf_plot(h(3),ni,'b');
irf_zoom(h(3),'y',[1 100]);
hold(h(3),'off');
set(h(3),'yscale','log');
ylabel(h(3),'n (cm^{-3})','Interpreter','tex');
irf_legend(h(3),{'n_{e}','n_{i}'},[0.5 0.9]);
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('eomni');
irf_spectrogram(h(4),specomni,'log');
irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(4),'E (eV)','fontsize',12,'Interpreter','tex');

h(5)=irf_panel('epadl');
irf_spectrogram(h(5),specpadl,'log');
irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_zoom(h(5),'y',[0 180]);
set(h(5),'ytick',[0 45 90 135 180]);
c_eval('ylabel(h(5),{''\theta (deg)'',''E > ? keV''},''fontsize'',12,''Interpreter'',''tex'');',Ethres/1000);

h(6)=irf_panel('epads');
irf_spectrogram(h(6),specpads,'log');
irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',12)
irf_zoom(h(6),'y',[0 180]);
set(h(6),'ytick',[0 45 90 135 180]);
c_eval('ylabel(h(6),{''\theta (deg)'',''E < ? keV''},''fontsize'',12,''Interpreter'',''tex'');',Ethres/1000);

c_eval('title(h(1),''MMS ?'');',ic);

irf_plot_axis_align(h(1:6));
irf_zoom(h(1:6),'x',tint);
set(h(1:6),'fontsize',12);

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%c_eval('print(''-dpng'',''-painters'',''-r600'',''MMS?overviewePADs.png'');',ic);