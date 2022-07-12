%% Load the data
tint = irf.tint('2021-11-03T12:10:00.000Z/2021-11-03T12:40:00.000Z');
B = solo.get_data('B_rtn_norm',tint);
scpot = solo.get_data('scpot', tint);
E = solo.get_data('e_rtn',tint);
Vi = solo.get_data('Vi_rtn',tint);
Ni = solo.get_data('Ni',tint);
Ne = solo.psp2ne(scpot);
Ti = solo.get_data('Ti',tint);
Eiflux = solo.get_data('pas_eflux',tint);
QF = solo.get_data('pas_qf',tint);

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

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

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
colormap(h(4),'jet'); ylabel(h(4),'f [Hz]')

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

irf_spectrogram(h(10),Eiflux); 
set(h(10), 'YScale', 'log'); 
ylabel(h(10),'E [eV]');
colormap(h(10),'jet');

set(h,'linewidth',1,'fontsize',15)
axis(h,'tight')
irf_zoom(h,'x',tint)
irf_plot_axis_align(h);

title(h(1),['SWA quality factor ' num2str(min(QF.data)) ' - ' num2str(max(QF.data))],'fontsize',15,'FontWeight','normal')
h(1).TitleHorizontalAlignment = 'right';