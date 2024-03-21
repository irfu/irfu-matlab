% Perform polarization analysis on burst mode electric and magnetic fields.
% Plots spectrograms, ellipticity, wave-normal angle, planarity, degree of
% polarization (DOP), phase speed, and normalized Poynting flux along B.
% Time selections should not be too long (less than 20 seconds), 
% otherwise the analysis will be very slow. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number

Tint = irf.tint('2017-06-19T09:43:19.00Z/2017-06-19T09:43:30.00Z');
Tint = irf.tint('2017-06-19T09:43:25.00Z/2017-06-19T09:43:28.00Z');
Tint = irf.tint('2017-07-06T00:54:10.00Z/2017-07-06T00:54:40.00Z');

%% Load data
Tintl = Tint+[-100 100];
R  = mms.get_data('R_gse',Tintl);
c_eval('Rxyz = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);

%% Polarization analysis
Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
B_SI=Bxyz.abs.data*1e-9;
Wce = e*B_SI/Me;
ecfreq = Wce/2/pi;
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = irf.ts_scalar(Bxyz.time,ecfreq);
ecfreq01 = irf.ts_scalar(Bxyz.time,ecfreq01);
ecfreq05 = irf.ts_scalar(Bxyz.time,ecfreq05);

polarization = irf_ebsp(Exyz,Bscm,Bxyz,Bxyz,Rxyz,[2 4000],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Bperp = polarization.bb_xxyyzzss(:,:,1)+polarization.bb_xxyyzzss(:,:,2);
Esum = polarization.ee_xxyyzzss(:,:,4);
Eperp = polarization.ee_xxyyzzss(:,:,1)+polarization.ee_xxyyzzss(:,:,2);
Epar = polarization.ee_xxyyzzss(:,:,3);
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);
% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;
vphperp = sqrt(Eperp./Bperp)*1e6;

% Remove points with very low B amplitutes
Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
pfluxz(removepts) = NaN;
vph(removepts) = NaN;
vphperp(removepts) = NaN;

%% Plot

h=irf_plot(8,'newfigure'); 
%h=irf_figure(540+ic,8);
xSize=750; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.88;
ywidth = 0.115;
set(h(1),'position',[0.08 0.975-ywidth xwidth ywidth]);
set(h(2),'position',[0.08 0.975-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.08 0.975-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.08 0.975-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.08 0.975-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.08 0.975-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.08 0.975-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.08 0.975-8*ywidth xwidth ywidth]);

h(1)=irf_panel('Bsum');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Bsum;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
irf_spectrogram(h(1),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(1),'(a)',[0.99 0.98],'color','w','fontsize',12)
hold(h(1),'on');
irf_plot(h(1),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(1),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(1),ecfreq01,'linewidth',1.5,'color','w')
hold(h(1),'off');
set(h(1),'yscale','log');
set(h(1),'ytick',[1e1 1e2 1e3]);
caxis(h(1),[-8 -1])
ylabel(h(1),'f (Hz)','fontsize',12);

h(2)=irf_panel('Esum');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Esum;
specrec.f_label='';
specrec.p_label={'log_{10}E^{2}','mV^2 m^{-2} Hz^{-1}'};
irf_spectrogram(h(2),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(2),'(b)',[0.99 0.98],'color','w','fontsize',12)
hold(h(2),'on');
irf_plot(h(2),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(2),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(2),ecfreq01,'linewidth',1.5,'color','w')
hold(h(2),'off');
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3]);
caxis(h(2),[-6 1])
ylabel(h(2),'f (Hz)','fontsize',12);

h(3)=irf_panel('ellipt');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=ellipticity;
specrec.f_label='';
specrec.p_label={'Ellipticity'};
irf_spectrogram(h(3),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(3),'(c)',[0.99 0.98],'color','w','fontsize',12)
hold(h(3),'on');
irf_plot(h(3),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(3),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(3),ecfreq01,'linewidth',1.5,'color','w')
hold(h(3),'off');
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3]);
caxis(h(3),[-1, 1])
ylabel(h(3),'f (Hz)','fontsize',12);

h(4)=irf_panel('thetak');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=thetak;
specrec.f_label='';
specrec.p_label={'\theta_{k}'};
irf_spectrogram(h(4),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
hold(h(4),'on');
irf_plot(h(4),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq01,'linewidth',1.5,'color','w')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[0, 90])
ylabel(h(4),'f (Hz)','fontsize',12);

h(5)=irf_panel('dop');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=dop;
specrec.f_label='';
specrec.p_label={'DOP'};
irf_spectrogram(h(5),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
hold(h(5),'on');
irf_plot(h(5),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq01,'linewidth',1.5,'color','w')
hold(h(5),'off');
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3]);
caxis(h(5),[0, 1])
ylabel(h(5),'f (Hz)','fontsize',12);

h(6)=irf_panel('planarity');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=planarity;
specrec.f_label='';
specrec.p_label={'planarity'};
irf_spectrogram(h(6),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',12)
hold(h(6),'on');
irf_plot(h(6),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(6),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(6),ecfreq01,'linewidth',1.5,'color','w')
hold(h(6),'off');
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3]);
caxis(h(6),[0, 1])
ylabel(h(6),'f (Hz)','fontsize',12);
  
h(7)=irf_panel('vph');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=vph;
specrec.f_label='';
specrec.p_label={'log_{10}E/B','m s^{-1}'};
irf_spectrogram(h(7),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','w','fontsize',12)
hold(h(7),'on');
irf_plot(h(7),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(7),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(7),ecfreq01,'linewidth',1.5,'color','w')
hold(h(7),'off');
set(h(7),'yscale','log');
set(h(7),'ytick',[1e1 1e2 1e3]);
caxis(h(7),[5 8])
ylabel(h(7),'f (Hz)','fontsize',12);

if 0
h(8)=irf_panel('poynting');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=pfluxz;
specrec.f_label='';
specrec.p_label={'S_{||}/|S|'};
irf_spectrogram(h(8),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','w','fontsize',12)
hold(h(8),'on');
irf_plot(h(8),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(8),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(8),ecfreq01,'linewidth',1.5,'color','w')
hold(h(8),'off');
set(h(8),'yscale','log');
set(h(8),'ytick',[1e1 1e2 1e3]);
caxis(h(8),[-1 1])
ylabel(h(8),'f (Hz)','fontsize',12);
end

hca=irf_panel('Epar');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Epar;
specrec.f_label='';
specrec.p_label={'log_{10}E_{||}^{2}','mV^2 m^{-2} Hz^{-1}'};
irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
hold(hca,'on');
irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-6 1])
ylabel(hca,'f (Hz)','fontsize',12);

% Remove grid and set background to grey
set(h(1:9),'xgrid','off','ygrid','off')
set(h(1:9),'Color',0.7*[1 1 1]);

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

colormap(h(1),'jet');
colormap(h(2),'jet');
colormap(h(3),bgrcmap);
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(7),'jet');
%colormap(h(8),bgrcmap);
%colormap(h(9),'jet');

c_eval('title(h(1),''MMS? Polarization Analysis'');',ic);
irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
set(h(1:8),'fontsize',12);
