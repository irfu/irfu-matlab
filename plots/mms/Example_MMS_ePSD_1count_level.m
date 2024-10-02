%%Example script for calculating the 1 count level for the electron
%distribution function measured by MMS.
%For more details on the various formulas used see the FPI instrument
%paper Pollock et al. 2016. In addition the Cluster RAPID calibration
%report is very useful as well, specifically Appendices A,B, and C
%( https://www2.mps.mpg.de/dokumente/projekte/cluster/rapid/caa_rap_cal.pdf )
%
%Written by Ahmad Lalti.
%% Define time interval and sampling rate of interest (a shock crossing)
clearvars
Tint = irf.tint('2023-04-24T03:49:00.000000000Z','2023-04-24T03:50:40.000000000Z');
Tu = irf.tint('2023-04-24T03:50:26.43Z','2023-04-24T03:50:35.02Z');
smplng = 'brst';
ic = 1;
%% Load 1count level from SDC
ePDist_ol = mms.db_get_ts(['mms1_fpi_' smplng '_l2_des-dist'],'mms1_des_avgf1counts_brst',Tint);
%Convert to SI units (because why not?!) - note that the convertto function
%doesn't work on this variable
un = split(ePDist_ol.siConversion,'>');
ePDist_ol = ePDist_ol*str2double(un{1});
ePDist_ol.units = un{2};
f1c_sdc = nanmean(ePDist_ol.data,1);
%% Sometimes it is not available, so estimate it yourself: f_1c = G*1, where
%G is the geometric factor.

%Load F_e and dF_e
ePDist = mms.get_data(['PDe_fpi_' smplng '_l2'],Tint,ic);
ePDistErr = mms.get_data(['PDERRe_fpi_' smplng '_l2'],Tint,ic);

ePDist = ePDist.convertto('s^3/m^6');
ePDistErr = ePDistErr.convertto('s^3/m^6');
%convert to counts and calculate the geometric factor (f = G*C)
C = ePDist;
C.data = ((ePDist.data./ePDistErr.data).^2);%from PSD to counts
G = ePDist.data./C.data;%Geometric factor of each instrument bin


%%my estimate for the 1 count level
ePDist_olm = ePDist;
ePDist_olm.data = ones(size(G)).*G;
fomni = ePDist_olm.omni;

%%By comparing my calculation of the 1 count level to that from the SDC for
%%various intervals in the solar wind I found that there is a factor of
%%around 1.2-1.4 difference between the two estimates, so I multiply
%%my estimate by 1.5 to be conservative (maybe try to calibrate yourself
%%for an interval in your region of interest? I did my
%%calculations in the solar wind)
wf = 1.5;
f1c_calc = nanmean(fomni.tlim(Tu).data,1);
f1c_calc_w = nanmean(fomni.tlim(Tu).data,1)*wf;
E_1c = nanmean(fomni.tlim(Tu).depend{1},1);

%% plot to compare the two 1count level

figure
loglog(E_1c,f1c_sdc,'k','linewidth',2,'Displayname','SDC')
hold on
loglog(E_1c,f1c_calc,'b','linewidth',2,'Displayname','Calc')
loglog(E_1c,f1c_calc_w,'r--','linewidth',2,'Displayname','Calc-weighted')
grid on
legend
ylabel('f_{1c} (s^3/m^6)')
xlabel('E (eV)')
set(gca,'FontSize',30)