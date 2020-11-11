% solo.correlate_probes_batch
%
% Script to get basic calibration coefficients for BIAS E and SCPOT
%
%% Load data
% One month, but we take +/- 3 days as a margin to better match with the
% nearby month
Tint = irf.tint('2020-10-28T00:00:00Z/2020-12-03T23:59:59Z'); 

VDC = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e', 'VDC', Tint);
EDC = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e', 'EDC', Tint);
QUAL = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e', 'QUALITY_FLAG', Tint);

% Remover thrusters, etc
iiBad = QUAL.data <2;
EDC.data(iiBad,:) = NaN; VDC.data(iiBad,:) = NaN; 

%% Step 1 - get d23
OutTS_step1 = solo.correlate_probes(VDC, EDC);
clf
irf_plot(OutTS_step1,'.-')
legend('d23 (Volts)', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
title('Step 1')

d23 = irf.ts_scalar(OutTS_step1.time,movmedian(OutTS_step1.data(:,1),11,'omitnan','Endpoints','fill'));

%% Step 2 - get k123
OutTS_step2 = solo.correlate_probes(VDC,EDC,d23);

irf_plot(OutTS_step2,'.-')
legend('d23 (Volts)', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
title('Step 2')

k123 = irf.ts_scalar(OutTS_step1.time,movmedian(OutTS_step2.data(:,3),15,'omitnan','Endpoints','fill'));

%% Step 3 - get d123
OutTS_l3 = solo.correlate_probes(VDC,EDC,d23,k123);

irf_plot(OutTS_l3,'.-')
legend('d23 (Volts)', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
title('Step 3')

d123 = irf.ts_scalar(OutTS_step1.time,movmedian(OutTS_l3.data(:,4),7,'omitnan','Endpoints','fill'));

%% Save result intoa monthly file
K123 = irf.ts_scalar(OutTS_step1.time,[k123.data, d123.data]);
save d23K123_20201110_november K123 d23

%% Validation plot 1
% downsample to 10 sec resolution to make plotting easier
Tstart = VDC.time.start+3600;
TstartS = Tstart.toUtc;
TstartS = [TstartS(1:11) '00:00:00Z'];
Tstart = EpochTT(TstartS);

Tstop = VDC.time.stop+3600;
TstopS = Tstop.toUtc;
TstopS = [TstopS(1:11) '00:00:00Z'];
Tstop = EpochTT(TstopS);

DT = 10; nSteps =  (Tstop-Tstart)/DT;
outTime = Tstart + ((1:nSteps)-0.5)*DT;

[DCE_SRF,PSP,SCPOT] = solo.vdccal(VDC,EDC);
DCE_SRF_10s = DCE_SRF.resample(outTime);
PSP_10s = PSP.resample(outTime);

h = irf_figure(948273,4,'reset');
irf_plot({DCE_SRF_10s,k123,d123,PSP_10s})

%% Validation plot 2 - reuires B from SOAR
B_RTN = solo.db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);
B_SRF = B_RTN ;
B_SRF.data(:,1:2) = -B_SRF.data(:,1:2);

E_sw = irf_e_vxb([-250 0 0],B_SRF);

h = irf_figure(948274,4,'reset');

hca = irf_panel('B');
hl = irf_plot(hca,B_SRF)
ylabel(hca,'B [nT]')

%hca = irf_panel('Ex');
%hl = irf_plot(hca,DCE_SRF.x);
%hold(hca,'on')
%hl = irf_plot(hca,E_sw.x,'r');
%hold(hca,'off')
%ylabel(hca,'Ex SRF [mV/m]')
%legend(hca,'E','-Vsw x B')

hca = irf_panel('Ey');
hl = irf_plot(hca,DCE_SRF_10s.y);
hold(hca,'on')
hl = irf_plot(hca,E_sw.y,'r');
hold(hca,'off')
ylabel(hca,'Ey SRF [mV/m]')

hca = irf_panel('Ez');
hl = irf_plot(hca,DCE_SRF_10s.z);
hold(hca,'on')
hl = irf_plot(hca,E_sw.z,'r');
hold(hca,'off')
ylabel(hca,'Ez SRF [mV/m]')

set(h(2:3), 'YLim', 4.99*[-1 1])


hca = irf_panel('PSP');
hl = irf_plot(hca,PSP_10s);
ylabel(hca,'PSP [V]')

irf_zoom(h,'x',Tint)

%% Load all month and combine
mm = {'january','february','march','april','may','june','july','august','september','october'};
K123 = []; d23 = [];
for m=3:10
  dd = load(['d23K123_20201110_' mm{m} '.mat']);
  tStart = EpochTT(sprintf('2020-%02d-01T00:00:00Z',m));
  tStop = EpochTT(sprintf('2020-%02d-01T00:00:00Z',m+1));
  if m==3, tStart = dd.d23.time(1); 
  elseif m==10, tStop = dd.d23.time(end); 
  end
  td23 = dd.d23.tlim(irf.tint(tStart,tStop));
  tK123 = dd.K123.tlim(irf.tint(tStart,tStop));
  if isempty(d23), d23 = td23; K123 = tK123; 
  else, d23 = d23.combine(td23); K123 = K123.combine(tK123);
  end
end

h = irf_figure(948278,1,'reset');
irf_plot(K123);
hold on
irf_plot(d23);
legend('k123, V1=k123*V23+del123','del123 (Volts)', 'd23 (Volts)')

%% BIAS change times
biasUpdateTimes = [...
  '2020-05-18T04:05:55Z';...
  '2020-06-14T00:00:00Z';...
  '2020-06-22T00:00:00Z';...
  '2020-07-06T00:00:00Z';...
  '2020-07-14T00:00:00Z';...
  '2020-07-20T07:39:41Z';...
  '2020-08-11T21:27:03Z';...
  '2020-09-05T00:00:00Z';...
  '2020-09-14T00:00:00Z'];

irf_plot(irf.ts_scalar(EpochTT(biasUpdateTimes),zeros(size(biasUpdateTimes,1),1)),'*')
legend('k123, V1=k123*V23+del123','del123 (Volts)', 'd23 (Volts)', 'Bias current update')

%% save the result
save d23K123_20201111 K123 d23