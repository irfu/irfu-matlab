% Plots E and B time series and of burst mode electric field in GSE
% coordinates and field-aligned coordinates. Plots spectrograms of parallel
% and perpendicular electric fields and fluctuating magnetic field.
% Written by D. B. Graham.

ic = 1; % Spacecraft number

Tint = irf.tint('2017-07-06T00:54:00.00Z/2017-07-06T00:54:45.00Z');
%Tint = irf.tint('2015-10-30T05:15:20.00Z/2015-10-30T05:16:20.00Z');
%Tint = irf.tint('2017-01-27T12:05:00.000Z/2017-01-27T12:06:00.000Z');


%% Load Data
c_eval('Bxyz=mms.get_data(''B_gse_brst_l2'',Tint,?);',ic);
c_eval('Exyz=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',ic);
c_eval('Bscm=mms.get_data(''B_gse_scm_brst_l2'',Tint,?);',ic);
c_eval('ne = mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ic);
magB = Bxyz.abs;
Bxyzmag = TSeries(Bxyz.time,[Bxyz.data magB.data]);

%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);

%% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm.time.epochUnix));
fmin = 20; fmax = 4100; %Hz
Exyzfachf = Exyzfac.filt(fmin,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,fmin,dfE,5);
Bscmfachf = Bscmfac.filt(fmin,0,dfB,5);

%% Wavelet transforms
nf = 100;
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[fmin fmax]);
Bwavelet = irf_wavelet(Bscm,'nf',nf,'f',[fmin fmax]);

%compress wavelet transform data 10 point average
nc = 100;
idx = nc/2:nc:length(Ewavelet.t)-nc/2;
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = 1:length(idx)
  Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specperpE=struct('t',Ewavelettimes);
specperpE.f=Ewavelet.f;
specperpE.p=Ewaveletx+Ewavelety;
specperpE.f_label='';
specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};

specparE=struct('t',Ewavelettimes);
specparE.f=Ewavelet.f;
specparE.p=Ewaveletz;
specparE.f_label='';
specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};


idx = nc/2:nc:length(Bwavelet.t)-nc/2;
Bwavelettimes = Bwavelet.t(idx);
Bwaveletx = zeros(length(idx),nf);
Bwavelety = zeros(length(idx),nf);
Bwaveletz = zeros(length(idx),nf);
for ii = 1:length(idx)
  Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specB=struct('t',Bwavelettimes);
specB.f=Bwavelet.f;
specB.p=Bwaveletx+Bwavelety+Bwaveletz;
specB.f_label='';
specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);

%% Plot Figure

h=irf_plot(7,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.86;
ywidth = 0.13;
set(h(1),'position',[0.10 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.10 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.10 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.10 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.10 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.10 0.97-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.10 0.97-7*ywidth xwidth ywidth]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyzmag);
ylabel(h(1),{'B (nT)'},'Interpreter','tex');
irf_zoom(h(1),'y',[-50 60]);
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}','|B|'},[0.98 0.12])
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',12)

h(2)=irf_panel('Elf');
irf_plot(h(2),Exyzfaclf);
ylabel(h(2),{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(2),{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.98 0.12])
irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',12)

h(3)=irf_panel('Ehf');
irf_plot(h(3),Exyzfachf);
ylabel(h(3),{'\delta E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(3),{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.98 0.12])
irf_legend(h(3),'(c)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(3),sprintf('f > %.1f Hz',fmin),[0.1 0.1],'color','k','fontsize',12)

h(4)=irf_panel('Especperp');
irf_spectrogram(h(4),specperpE,'log');
hold(h(4),'on');
irf_plot(h(4),Flh,'color','k','LineWidth',1.5)
irf_plot(h(4),Fce,'color','r','LineWidth',1.5)
irf_plot(h(4),Fpp,'color','b','LineWidth',1.5)
hold(h(4),'off');
irf_legend(h(4),'(d)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(h(4),'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(h(4),'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(h(4),'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(h(4),[-6 2]);
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(4),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(5)=irf_panel('Especpar');
irf_spectrogram(h(5),specparE,'log');
hold(h(5),'on');
irf_plot(h(5),Flh,'color','k','LineWidth',1.5)
irf_plot(h(5),Fce,'color','r','LineWidth',1.5)
irf_plot(h(5),Fpp,'color','b','LineWidth',1.5)
hold(h(5),'off');
irf_legend(h(5),'(e)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(h(5),'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(h(5),'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(h(5),'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(h(5),[-6 2]);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(6)=irf_panel('Bscmhf');
irf_plot(h(6),Bscmfachf);
ylabel(h(6),{'\delta B (nT)'},'Interpreter','tex');
irf_legend(h(6),{'B_{\perp 1}','B_{\perp 2}','B_{||}'},[0.98 0.12])
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(6),sprintf('f > %.1f Hz',fmin),[0.1 0.1],'color','k','fontsize',12)

h(7)=irf_panel('Bspec');
irf_spectrogram(h(7),specB,'log');
hold(h(7),'on');
irf_plot(h(7),Flh,'color','k','LineWidth',1.5)
irf_plot(h(7),Fce,'color','r','LineWidth',1.5)
irf_plot(h(7),Fpp,'color','b','LineWidth',1.5)
hold(h(7),'off');
irf_legend(h(7),'(g)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(h(7),'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(h(7),'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(h(7),'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(h(7),[-8 -1]);
set(h(7),'yscale','log');
set(h(7),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(7),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

%load('caa/cmap.mat');
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(7),'jet');

c_eval('title(h(1),''MMS?'')',ic);

irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',Tint);
set(h(1:7),'fontsize',12);