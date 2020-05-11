function [iPhoto] = pl_sweep_batch()
%SOLO.PL_SWEEP_BATCH  Plot all BIAS sweeps

YY = 2020; MM = 2:4;

iPhoto = [];
for im=1:length(MM)
  dataDir = sprintf('/Volumes/solarorbiter/remote/data/BIA/%d/%02d/',YY,MM(im));
  d = dir([dataDir 'solo_L1_rpw-bia-sweep-*.cdf']);
  for i=1:length(d)
    iPhoto_out = solo.pl_sweep([d(i).folder '/' d(i).name],1);
    if isempty(iPhoto), iPhoto = iPhoto_out;
    else, iPhoto = iPhoto.combine(iPhoto_out);
    end
  end
end

%% Plot
h = irf_figure(93766,1,'reset');
irf_plot(h,iPhoto/1000,'.-')
Tint = irf.tint('2020-02-01T00:00:00Z/2020-05-01T00:00:00Z');
irf_zoom(h,'x',Tint)

ylabel(h,'I saturation [uA]')

legend(h,'V1','V2','V3','Location','southwest')