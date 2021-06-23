function [iPhoto] = pl_sweep_batch()
%SOLO.PL_SWEEP_BATCH  Plot all BIAS sweeps

YY = 2020; MM = 6;

iPhoto = [];
for im=1:length(MM)
  dataDir = sprintf('/Volumes/solarorbiter/remote/data/BIA/%d/%02d/',YY,MM(im));
  d = dir([dataDir 'solo_L1_rpw-bia-sweep_*.cdf']);
  for i=1:length(d)
    %if MM(im)==5 && i<15, continue, end
    iPhoto_out = solo.pl_sweep([d(i).folder '/' d(i).name],1);
    if isempty(iPhoto), iPhoto = iPhoto_out;
    else, iPhoto = iPhoto.combine(iPhoto_out);
    end
  end
end

%% Plot
h = irf_figure(93766,1,'reset');
irf_plot(h,iPhoto_new/1000,'.-')
Tint = irf.tint(sprintf('2020-%02d-01T00:00:00Z/2020-%02d-01T00:00:00Z',MM(1),MM(end)+1));
irf_zoom(h,'x',Tint)

ylabel(h,'I saturation [uA]')

legend(h,'V1','V2','V3','Location','southwest')