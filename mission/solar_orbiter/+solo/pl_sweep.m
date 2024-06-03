function [iPhoto_out] = pl_sweep(fName, hardcopyFlag)
%SOLO.PL_SWEEP  Plots a BIAS sweep
%
% [iPhoto] = solo.pl_sweep(fileName, [hardcopyFlag])
%
% Input:
%     hardCopyFlag - Logical. Whether to produce a plot file in the current
%                    directory. Default=false.
%
% Output:
%     iPhoto - photo-saturation current
%
% NOTE: This function is called by cron jobs for batch processing BIAS sweep
% plots on brain/spis at IRFU.

iPhoto = NaN(1,3);

PHOTO_LIM = [-55 -5];

if nargin<2, hardcopyFlag = false; end
dObj = dataobj(fName);
[~,fName,~] = fileparts(fName);

V=get_ts(dObj,'V'); time = V.time(1);
I=get_ts(dObj,'BIAS_SWEEP_CURRENT');
AntFlag = get_ts(dObj,'ANT_FLAG');

h = irf_figure(93485,2,'reset');
h = irf_plot({I*1e-3,V});
ylabel(h(1),'I [\mu A]')
title(h(1),fName,'Interpreter','none')

if hardcopyFlag

  set(gcf,'InvertHardcopy','off','PaperPositionMode','auto')
  print(gcf,'-dpng','-painters','-r150',[fName '.png'])
end

%%
h = irf_figure(93486,1,'reset');
for iProbe=1:3
  ii = AntFlag.data==iProbe;
  dataI = I.data(ii,iProbe);
  dataV = V.data(ii,iProbe);
  if isempty(dataI) & isempty(dataV)
    % NOTE: A sweep dataset may lack data for one probe, which then results in
    %       empty data arrays which can not be truncated.
    % Example: solo_L1_rpw-bia-sweep-cdag_20240522T020342-20240529T020451_V02.cdf
    irf.log('w', sprintf('"%s" does not contain any data for antenna %i.', fName, iProbe))
  else
    %Remove initial jumps
    dataI(1:2) = []; dataV(1:2) = [];
  end

  %Determine photo saturation current
  iiPhoto = dataV>PHOTO_LIM(1) & dataV<PHOTO_LIM(2);
  if any(iiPhoto)
    iPhoto(iProbe) = median(dataI(iiPhoto));
    stdPhoto = std(dataI(iiPhoto));
    if (stdPhoto > 5*1000) || (iPhoto(iProbe) > 0), iPhoto(iProbe) = NaN; end
  else, stdPhoto = NaN;
  end

  plot(h, dataV, dataI/1000,'.-')
  irf_legend(sprintf('I_{sat}%d = %.1f uA (std = %.1f)',iProbe,iPhoto(iProbe)/1000, stdPhoto/1000) ,[0.8, 0.4-0.1*iProbe]);
  switch iProbe
    case 1, hold(h,'on')
    case 3, hold(h,'off')
  end
end
ylabel(h,'I [\mu A]')
xlabel(h,'U [V]')
grid(h,'on')
set(h,'XLim',[-70 70])
legend(h,'V1','V2','V3','Location','northwest')
title(h(1),fName,'Interpreter','none')


iPhoto_out = irf.ts_scalar(time,iPhoto);
iPhoto_out.units = 'nA';
iPhoto_out.name = 'I saturation';

%%
%{
for ip=1:3
  %sprintf('Probe %d',ip)
  ii = find(I.data(:,ip) == 0);
  iJump = find(diff(ii)>1);
  iStart = ii(1);
  for ij = 1:length(iJump)+1
    if ij > length(iJump), iStop = ii(end);
    else, iStop = ii(iJump(ij));
    end
    if length(I.data(iStart:iStop,ip))>128
      %sprintf('blanking %d-%d',iStart,iStop)
      I.data(iStart:iStop,ip) = NaN;
      V.data(iStart:iStop,ip) = NaN;
    end
    if ij > length(iJump), continue, end
    iStart = ii(iJump(ij)+1);
  end
end


%%
h = irf_figure(93486,1,'reset');
plot(h,V.data, I.data/1000,'.-')
ylabel(h,'I [\mu A]')
xlabel(h,'U [V]')
grid(h,'on')
set(h,'XLim',[-70 70])
legend(h,'V1','V2','V3','Location','northwest')
title(h(1),fName,'Interpreter','none')
%}


if hardcopyFlag

  set(gcf,'InvertHardcopy','off','PaperPositionMode','auto')
  print(gcf,'-dpng','-painters','-r150',[fName '_IV.png'])
end