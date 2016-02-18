% Scipt to plot spin rate

tStart = '2015-04-04T05:00:00.000000Z';
dtSec = 4*3600;
defAttFile = '/Users/yuri/Data/mms/ancillary/mms1/defatt/MMS1_DEFATT_2015093_2015094.V00';

%% Load DEFATT and computer phase at 1 sec cadence
defatt = mms_load_ancillary(defAttFile,'defatt');
epoch0 = EpochTT(tStart);
t = int64((1:(dtSec))'*1e9) + epoch0.epoch;
pha = mms_defatt_phase(defatt,t);
phaUnwrapped = unwrap(pha.data*pi/180);

%% Plot
figure
irf_plot(irf.ts_scalar(EpochTT(t(2:end)),diff(phaUnwrapped)/2/pi*60))
ylabel('Spin Rate [rpm]')