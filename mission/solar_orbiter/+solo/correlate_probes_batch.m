DPATH = '/Volumes/solo/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/latest/';
d = dir([DPATH '*.cdf']);

probeCorr = [];
for i = 82:length(d) % 7 - March, 82 - June 1
  sprintf('Processing %s', d(i).name)
  tt = solo.correlate_probes([d(i).folder filesep d(i).name]);
  if isempty(probeCorr), probeCorr = tt;
  else, probeCorr = probeCorr.combine(tt);
  end
end

%%
irf_plot(probeCorr)
probeCorr.data(probeCorr.data(:,1)>1.3,:) = NaN;
probeCorr.data(probeCorr.data(:,1)<1.1,:) = NaN;
probeCorr.data(probeCorr.data(:,3)<0.5,:) = NaN;
irf_plot(probeCorr,'o-')
legend('k1, V1=k1*V3+d13','d13 (Volts)','k2, V2=k2*V3+d23','d23 (Volts)')

%%

d23 = irf.ts_scalar(probeCorr.time,movmedian(probeCorr.data(:,5),11,'omitnan','Endpoints','fill'));
idx = find(~isnan(d23.data));
d23 = irf.ts_scalar(d23.time(idx), d23.data(idx));
%%
DPATH = '/Volumes/solo/data_irfu/SOLO_L2_LFR-SURV-CWF-E_without_sweeps/latest/';
d = dir([DPATH '*.cdf']);

probeCorrNew = [];
for i = 82:length(d)
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

%% Save update data
K123_new = K123;
d23_new = d23;
load d23K123
K123 = K123_new.combine(K123); d23 = d23_new.combine(d23);
save d23K123 K123 d23