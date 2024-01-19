%% Load the data
Tint = irf.tint('2021-08-22T02:00:00.00Z/2021-08-22T02:40:00.00Z');
B = solo.get_data('B_rtn_norm',Tint);
scpot = solo.get_data('scpot', Tint);
E = solo.get_data('e_rtn',Tint);
Vi = solo.get_data('Vi_rtn',Tint);
Ni = solo.get_data('Ni',Tint);
Ne = solo.psp2ne(scpot);
%Ti = solo.get_data('Ti',tint); % Scalar temperature
Ti = solo.get_data('Ti_fac',Tint);
Eiflux = solo.get_data('pas_eflux',Tint);
QF = solo.get_data('pas_qf',Tint);

%% Calculate ellipticity and wavelet Bsum
b0 = irf_filt(B,0,0.01,[],3);
fs_B = 1 / (B.time(2)-B.time(1));
ebsp2 =irf_ebsp([],B,[],b0,[],[0.05 0.5*fs_B],...
  'polarization','fac','fullB=dB');
frequency = ebsp2.f;
time = ebsp2.t;
Bsum = ebsp2.bb_xxyyzzss(:,:,4);
ellipticity = ebsp2.ellipticity;
dop = ebsp2.dop;
dopthresh = 0.7;
removepts = find(dop < dopthresh);
ellipticity(removepts) = NaN;
msk = ellipticity;
msk(~isnan(msk)) = 1;
msk(isnan(msk)) = 0;
msk_denoise = bwareaopen(msk,8);
ellipticity(msk_denoise==0) = NaN;
spec_ellipticity=struct('t',time);
spec_ellipticity.f=frequency;
spec_ellipticity.p=ellipticity;
spec_ellipticity.f_label='f [Hz]';
spec_ellipticity.p_label={'Elipticity'};

specrec_Bsum=struct('t',time);
specrec_Bsum.f=frequency;
specrec_Bsum.p=Bsum;
specrec_Bsum.f_label='';
specrec_Bsum.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};

%% Plot the data
bgrcmap=irf_colormap('bluered');
spacmap=irf_colormap('space');

h = irf_figure(10);
irf_plot(h(1),B.abs);
ylabel(h(1),'|B| [nT]','Interpreter','tex')

irf_plot(h(2),B)
ylabel(h(2),'B_{rtn} [nT]','Interpreter','tex')
irf_legend(h(2),{'B_r';'B_t';'B_n'},[1.02 0.98],'fontsize',15)

E_hp=irf_filt(E,0.5,0,[],3);
irf_plot(h(3),E_hp)
ylabel(h(3),'E_{rtn} [nT]','Interpreter','tex')
irf_legend(h(3),{'E_r';'E_t';'E_n'},[1.02 0.98],'fontsize',15)

irf_spectrogram(h(4),specrec_Bsum);
set(h(4), 'YScale', 'log')
colormap(h(4),spacmap); ylabel(h(4),'f [Hz]')

irf_spectrogram(h(5),spec_ellipticity,'log');
caxis(h(5),[-1 1])
set(h(5), 'YScale', 'log');
colormap(h(5),bgrcmap)

irf_plot(h(6),Ni)
hold(h(6))
irf_plot(h(6),Ne,'color','b')
ylabel(h(6),'N [cc]','Interpreter','tex')
irf_legend(h(6),{'N PAS';'N RPW'},[1.02 0.98],'fontsize',15)

irf_plot(h(7),Vi.x,'color','k')
ylabel(h(7),'V_{r} [km/s]','Interpreter','tex')

irf_plot(h(8),Vi.y,'color','b')
hold(h(8))
irf_plot(h(8),Vi.z,'color','r')
ylabel(h(8),'V_{tn} [km/s]','Interpreter','tex')
irf_legend(h(8),{'';'V_t';'V_n'},[1.02 0.98],'fontsize',15)

irf_plot(h(9),Ti)
ylabel(h(9),'T_{i} [eV]','Interpreter','tex')
if size(Ti.data,2)>1
  irf_legend(h(9),{'T_{||}';'T_{\perp 1}';'T_{\perp 2}'},[1.02 0.98],'fontsize',15)
end

irf_spectrogram(h(10),Eiflux);
set(h(10), 'YScale', 'log');
ylabel(h(10),'E [eV]');
colormap(h(10),spacmap);


set(h,'linewidth',1,'fontsize',15)
axis(h,'tight')
irf_zoom(h,'x',Tint)
irf_plot_axis_align(h);

set(h(10), 'YLim', [300 5000]);

title(h(1),['SWA quality factor ' num2str(min(QF.data)) ' - ' num2str(max(QF.data))],'fontsize',15,'FontWeight','normal')
h(1).TitleHorizontalAlignment = 'right';