function [iPhoto] = pl_sweep_batch(YY,MM)
%SOLO.PL_SWEEP_BATCH  Plot all BIAS sweeps for one month
%
% [iPhoto] = pl_sweep_batch(YYYY,MM)

%YY = 2023; MM = 3;

iPhoto = [];
for im=1:length(MM)
  dataDir = sprintf('/Volumes/solo/remote/data/BIA/%d/%02d/',YY,MM(im));
  d = dir([dataDir 'solo_L1_rpw-bia-sweep*.cdf']);
  for i=1:length(d)
    %if MM(im)==5 && i<15, continue, end
    iPhoto_out = solo.pl_sweep([d(i).folder '/' d(i).name],1);
    if all(isnan(iPhoto_out.data)), continue, end
    if isempty(iPhoto), iPhoto = iPhoto_out;
    else, iPhoto = iPhoto.combine(iPhoto_out);
    end
  end
end

%% Plot
h = irf_figure(93766,1,'reset');
irf_plot(h,iPhoto/1000,'.-')
YY2 = YY; MM2 = MM(end)+1;
if MM2>12, YY2 = YY2+1; MM2 = MM2-12; end
Tint = irf.tint(sprintf('%d-%02d-01T00:00:00Z/%d-%02d-01T00:00:00Z',YY,MM(1),YY2,MM2));
irf_zoom(h,'x',Tint)

ylabel(h,'I saturation [uA]')

legend(h,'V1','V2','V3','Location','southwest')

hardcopyFlag = true;
if hardcopyFlag
  fName = ['solo_Iph_sat_' irf_fname(Tint,5)];
  set(gcf,'InvertHardcopy','off','PaperPositionMode','auto')
  print(gcf,'-dpng','-painters','-r150',[fName '.png'])
end