% Plots E time series and of burst mode electric field in spacecraft DSL
% coordinates and field-aligned coordinates. Plots spectrograms of parallel
% and perpendicular electric fields. 
% Written by D. B. Graham.

% Select spacecraft number: 1--4
ic = 3;

tint = irf.tint('2015-08-28T12:53:15.00Z/2015-08-28T12:54:40.00Z');
c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

%load brst mode interval locally
%tmpDataObj = dataobj('data/mms3_edp_brst_ql_dce2d_20150828125314_v0.2.0.cdf');
%Exyz = get_variable(tmpDataObj,'mms3_edp_dce_xyz_dsl');
%Exyz = mms.variable2ts(Exyz);
%tmpDataObj = dataobj('data/mms3_dfg_srvy_ql_20150828_v0.0.3.cdf');
%Bxyz = get_variable(tmpDataObj,'mms3_dfg_srvy_dmpa');
%Bxyz = mms.variable2ts(Bxyz);

%Calculate electron gyrofrequencies
ecfreq = (1.6e-19)*Bxyz.abs.data*1e-9/(9.1e-31*2*pi);
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = TSeries(Bxyz.time,ecfreq);
ecfreq01 = TSeries(Bxyz.time,ecfreq01);
ecfreq05 = TSeries(Bxyz.time,ecfreq05);

Bxyz = Bxyz.tlim(tint);
Exyz = Exyz.tlim(tint);
Bxyz = Bxyz.resample(Exyz);

% Perform coordinate transformation into Field-Aligned Coordinates (FAC)
SCpos = [0 1 0];

Bmag = Bxyz.abs.data;
Rpar = Bxyz.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos);
Rmag   = irf_abs(Rperpy,1);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag   = irf_abs(Rperpx,1);
Rperpx = Rperpx./[Rmag Rmag Rmag];

Epar = dot(Rpar,Exyz.data,2);
Eperp = dot(Rperpx,Exyz.data,2);
Eperp2 = dot(Rperpy,Exyz.data,2);

Efac = TSeries(Exyz.time,[Eperp Eperp2 Epar],'to',1);

Ewavelet = irf_wavelet(Efac,'returnpower',1,'cutedge',1,'nf',40);


% Begin plotting
h=irf_plot(5,'newfigure'); 

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Exyz');
irf_plot(h(2),Exyz);
ylabel(h(2),'E_{DSL} (mV m^{-1})','Interpreter','tex');
irf_legend(h(2),{'E_{x}','E_{y}','E_{z}'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Efac');
irf_plot(h(3),Efac);
ylabel(h(3),'E_{FAC} (mV m^{-1})','Interpreter','tex');
irf_legend(h(3),{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.5 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('Eparwavelet');
  specrec=struct('t',Ewavelet.t);
    specrec.f=Ewavelet.f;
    specrec.p=Ewavelet.p{1,3}; % Power of parallel electric field
    specrec.f_label='';
    specrec.p_label={'E_{||}^2','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(4),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
%overplot electron gyrofrequencies, fce, 0.5fce, and 0.1fce
hold(h(4),'on');
irf_plot(h(4),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq01,'linewidth',1.5,'color','w')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[-6 0])
ylabel(h(4),'f (Hz)','fontsize',12);

h(5)=irf_panel('Eperpwavelet');
  specrec=struct('t',double(Ewavelet.t));
    specrec.f=Ewavelet.f;
    specrec.p=Ewavelet.p{1,1}+Ewavelet.p{1,2}; % Power of perpendicular electric fields
    specrec.f_label='';
    specrec.p_label={'E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(5),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
%overplot electron gyrofrequencies, fce, 0.5fce, and 0.1fce
hold(h(5),'on');
irf_plot(h(5),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq01,'linewidth',1.5,'color','w')
hold(h(5),'off');
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3]);
caxis(h(5),[-6 0])
ylabel(h(5),'f (Hz)','fontsize',12);

irf_zoom(h(1:5),'x',tint);
set(h(1:5),'fontsize',12);
irf_plot_axis_align(h(1:5));

pictitle = strcat('MMS',num2str(ic),' - Burst mode electric field data');
title(h(1),pictitle);

set(gcf,'paperpositionmode','auto')
print('-dpng','-painters','-r500','burstmodeEfields.png');