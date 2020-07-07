DPATH = '/Volumes/solarorbiter/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/2020-06-17_10.44/';
d = dir([DPATH '*.cdf']);
DPATH = '/Volumes/solarorbiter/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/2020-06-24_11.51/';
d1 = dir([DPATH '*.cdf']);
d = [d; d1];
probeCorr = [];
for i = 6:length(d)
  sprintf('Processing %s', d(i).name)
  tt = solo.correlate_probes([d(i).folder filesep d(i).name]);
  if isempty(probeCorr), probeCorr = tt;
  else, probeCorr = probeCorr.combine(tt);
  end
end

%%
irf_plot(probeCorr)
probeCorr.data(probeCorr.data(:,1)>1.8,:) = NaN;
probeCorr.data(probeCorr.data(:,1)<1.1,:) = NaN;
probeCorr.data(probeCorr.data(:,3)<0.5,:) = NaN;
irf_plot(probeCorr,'o-')
legend('k1, V1=k1*V3+d13','d13 (Volts)','k2, V2=k2*V3+d23','d23 (Volts)')

%%

d23 = irf.ts_scalar(probeCorr.time,movmedian(probeCorr.data(:,5),11,'omitnan','Endpoints','fill'));
idx = find(~isnan(d23.data));
d23 = irf.ts_scalar(d23.time(idx), d23.data(idx));
%%
DPATH = '/Volumes/solarorbiter/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/2020-06-17_10.44/';
d = dir([DPATH '*.cdf']);
DPATH = '/Volumes/solarorbiter/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/2020-06-24_11.51/';
d1 = dir([DPATH '*.cdf']);
d = [d; d1];
probeCorrNew = [];
for i = 6:length(d)
  sprintf('Processing %s', d(i).name)
  tt = solo.correlate_probes([d(i).folder filesep d(i).name],d23);
  if isempty(probeCorrNew), probeCorrNew = tt;
  else, probeCorrNew = probeCorrNew.combine(tt);
  end
end

%%

probeCorrNew.data(probeCorrNew.data(:,1)<0.9,:) = NaN;

K123 = irf.ts_scalar(probeCorrNew.time,movmedian(probeCorrNew.data,5,'omitnan','Endpoints','fill'));
idx = find(~isnan(K123.data(:,1)));
K123 = irf.ts_scalar(K123.time(idx), K123.data(idx,:));

%%
save d23K123 K123 d23