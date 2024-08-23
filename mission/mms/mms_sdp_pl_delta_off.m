function res = mms_sdp_pl_delta_off(Tint,mmsId)
%MMS_SDP_PL_DELTA_OFF  Delta offsets and FPI comparison
%
% res = MMS_SDP_PL_DELTA_OFF(Tint,mmsId)

%%
fPre = sprintf('mms%d_edp_fast_l2a_dce2d',mmsId);
Bmask = mms.db_get_ts(fPre,...
  sprintf('mms%d_edp_bitmask_fast_l2a',mmsId), Tint);
if isempty(Bmask), res = []; return, end
bitmask = Bmask.data;
epochFull = Bmask.time.ttns;
Es12 = mms.db_get_ts(fPre,...
  sprintf('mms%d_edp_espin_p12_fast_l2a',mmsId), Tint);
Es34 = mms.db_get_ts(fPre,...
  sprintf('mms%d_edp_espin_p34_fast_l2a',mmsId), Tint);
epochSpin = Es34.time.ttns;

mskAsponOn = bitand(bitmask(:,2),64) > 0;
mskAsponOn = mskAsponOn | (bitand(bitmask(:,3),64) > 0);

% XXX hack
if ~any(mskAsponOn)
  TAspOn = Tint.stop+(-3*3600);
  mskAsponOn = epochFull>TAspOn.ttns;
end

% Resample bitmask
spinSize = int64(20000000000);
mskAsponOnSpin = zeros(size(epochSpin));
ints = find_on(mskAsponOn);
for i=1:size(ints,1)
  mskAsponOnSpin((epochSpin> epochFull(ints(1,1))-spinSize/2) & ...
    (epochSpin<= epochFull(ints(1,2))+spinSize/2)) = 1;
end
mskAsponOnSpin=logical(mskAsponOnSpin);

%
Es12AspocOff =  Es12(~mskAsponOnSpin); Es12AspocOn =  Es12(mskAsponOnSpin);
Es34AspocOff =  Es34(~mskAsponOnSpin); Es34AspocOn =  Es34(mskAsponOnSpin);
DeltaAspocOff = Es12AspocOff - Es34AspocOff;
DeltaAspocOn = Es12AspocOn - Es34AspocOn;

%% FPI
B = mms.get_data('B_dmpa_srvy',Tint,mmsId);
if isempty(B), B = mms.get_data('B_dmpa_dfg_srvy_ql',Tint,mmsId); end

% First try L2
Vifpi = mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,mmsId);
if ~isempty(Vifpi)
  Nifpi = mms.get_data('Ni_fpi_fast_l2',Tint,mmsId);
  Tifpi = mms.get_data('Tsi_fpi_fast_l2',Tint,mmsId);
else % QL
  Vifpi = mms.get_data('Vi_gse_fpi_ql',Tint,mmsId);
  if ~isempty(Vifpi)
    Nifpi = mms.get_data('Ni_fpi_ql',Tint,mmsId);
    Tifpi = mms.get_data('Tsi_fpi_ql',Tint,mmsId);
  else %isempty(Vifpi)
    res = []; return
    Nifpi = []; Tifpi = []; Efpi = []; IonSpc = []; %#ok<UNRCH>
  end
end

if ~isempty(Vifpi)
  Efpi = irf_e_vxb(Vifpi,B.resample(Vifpi));

  mmsIdS = sprintf('mms%d',mmsId);
  fPref = [mmsIdS '_fpi_fast_ql_dis'];
  iEnSp_pX = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_pX'], Tint); %positiveX
  if isempty(iEnSp_pX), IonSpc = [];
  else
    iEnSp_pY = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_pY'], Tint); %positiveY
    iEnSp_pZ = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_pZ'], Tint); %positiveZ
    iEnSp_mX = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_mX'], Tint); %positiveX
    iEnSp_mY = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_mY'], Tint); %positiveY
    iEnSp_mZ = mms.db_get_variable(fPref,[mmsIdS '_dis_energySpectr_mZ'], Tint); %positiveZ

    [~,energy] = hist([log10(10),log10(30e3)],32); energy = 10.^energy;
    IonSpc = struct('t',irf_time(iEnSp_pX.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    IonSpc.f = energy; % 0:31 - energy levels
    IonSpc.p = iEnSp_pX.data+iEnSp_pY.data+iEnSp_pZ.data+...
      iEnSp_mX.data+iEnSp_mY.data+iEnSp_mZ.data;%data matrix
    IonSpc.f_label = 'E [eV]';
    IonSpc.p_label = 'au';
    IonSpc.plot_type = 'log';
  end

  VHplushpca = mms.get_data('Vhplus_dbcs_hpca_srvy_l2',Tint,mmsId);
  if isempty(VHplushpca), Ehpca = [];
  else, Ehpca = irf_e_vxb(VHplushpca,B.resample(VHplushpca));
  end
end
PSP = mms.db_get_ts(sprintf('mms%d_edp_fast_l2_scpot',mmsId),...
  sprintf('mms%d_edp_psp_fast_l2',mmsId),Tint);

%% resample to 1 min
epoch1min = fix(Tint.start.epochUnix/60)*60:20:ceil(Tint.stop.epochUnix/60)*60;
%epoch1min = epoch1min*60;
Epoch1min = EpochUnix(epoch1min);

if ~isempty(PSP)
  PSPR = PSP.resample(Epoch1min,'median');
end
if ~isempty(Es12AspocOff), Es12AspocOffR = Es12AspocOff.resample(Epoch1min,'median');
else, Es12AspocOffR = Es12AspocOff;
end
if ~isempty(Es12AspocOn), Es12AspocOnR = Es12AspocOn.resample(Epoch1min,'median');
else, Es12AspocOnR = Es12AspocOn;
end
if ~isempty(Es34AspocOff), Es34AspocOffR = Es34AspocOff.resample(Epoch1min,'median');
else, Es34AspocOffR = Es34AspocOff;
end
if ~isempty(Es34AspocOn), Es34AspocOnR = Es34AspocOn.resample(Epoch1min,'median');
else, Es34AspocOnR = Es34AspocOn;
end
if ~isempty(Efpi)
  EfpiR = Efpi.resample(Epoch1min,'median');
end
if isempty(DeltaAspocOff), DeltaAspocOffR = DeltaAspocOff;
else, DeltaAspocOffR = DeltaAspocOff.resample(Epoch1min,'median');
end
if isempty(DeltaAspocOn), DeltaAspocOnR = DeltaAspocOn;
else, DeltaAspocOnR = DeltaAspocOn.resample(Epoch1min,'median');
end

idxMSH = [];
if  ~isempty(Nifpi)
  NifpiR = Nifpi.resample(Epoch1min,'median');
  VifpiR = Vifpi.resample(Epoch1min,'median');
  idxMSH = NifpiR.data>5 & VifpiR.x.data>-200;
  [~,idxTmp,idxTmp2] = intersect(DeltaAspocOffR.time.epoch,NifpiR.time.epoch(idxMSH));
  DeltaAspocOffRes = DeltaAspocOffR(idxTmp);
  % XXX check if NifpiR exists
  NifpiRes = NifpiR(idxTmp2);
  PSPRes = PSPR(idxTmp2);
end

if isempty(Ehpca), EhpcaR = [];
else, EhpcaR = Ehpca.resample(Epoch1min,'median');
end

%% Raw data figure
myCols = [[0 0 0];[.3 .3 .3];[0 0 1];[0.3010 0.7450 0.9330];[1 0 1];[1 0 0]];
if 1
  h = irf_figure(93,9,'reset');
  cmap = irf_colormap('space'); colormap(cmap)

  hca = irf_panel('B');
  irf_plot(hca,B);
  if any(B.abs.data>100), set(hca,'YLim',[-99 99]), end
  ylabel(hca,'B [nT]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

  hca = irf_panel('DIS-spectrogram');
  if ~isempty(IonSpc)
    irf_spectrogram(hca,IonSpc,'log','donotfitcolorbarlabel');
  end
  set(hca,'YScale','log')
  hold(hca,'on')
  irf_plot(hca,Tifpi)
  hold(hca,'off')
  ylabel(hca,'E(i) [eV]'), set(hca,'YTick',[10 100 1000 10000])

  irf_plot_axis_align(h);

  hca = irf_panel('Vfpi');
  irf_plot(hca,Vifpi);
  ylabel(hca,'Vi [km/s]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

  hca = irf_panel('Ni');
  irf_plot(hca,Nifpi);
  ylabel(hca,'Ni [cc]')
  % ScPot
  hca = axes('Position',get(hca,'Position'));
  hl = irf_plot(hca,PSP); set(hl,'Color','b')
  ylabel(hca,'PSP [V]')
  set(hca, 'YAxisLocation','right');
  set(hca,'Color','none'); % color of axis
  set(hca,'XColor','b','YColor','b');
  irf_zoom(hca,'x',Tint), xlabel(hca,'')
  set(hca,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required

  hca = irf_panel('Delta_off'); set(hca,'ColorOrder',myCols)
  irf_plot(hca,{DeltaAspocOffR.x,DeltaAspocOnR.x,DeltaAspocOffR.y,DeltaAspocOnR.y},'comp')

  hold(hca,'on')
  irf_plot(hca,DeltaAspocOffRes,'.')
  hold(hca,'off')
  ylabel(hca,'Delta [mV/m]')
  set(hca,'YLim',[-1.4 1.4])
  set(hca,'ColorOrder',[[0 0 0];[0 0 1]]);
  irf_legend(hca,{'X','Y'},[0.95, 0.95])

  leg = {'12','12ao','34','34ao','dis'};
  hca = irf_panel('Ex'); set(hca,'ColorOrder',myCols)
  irf_plot(hca,{Es12AspocOffR.x,Es12AspocOnR.x,...
    Es34AspocOffR.x,Es34AspocOnR.x},'comp');
  if ~isempty(Efpi)
    hold(hca,'on'), irf_plot(hca,EfpiR.x,'m.'), hold(hca,'off')
  end
  if ~isempty(Ehpca)
    hold(hca,'on'), irf_plot(hca,EhpcaR.x,'r.'), hold(hca,'off'), leg = [leg 'hpca'];
  end
  ylabel(hca,'Ex [mV/m]'), irf_legend(hca,leg,[0.95, 0.95])

  hca = irf_panel('Fpi_x'); set(hca,'ColorOrder',myCols)
  if ~isempty(Efpi)
    irf_plot(hca,{EfpiR.x-Es12AspocOffR.x,EfpiR.x-Es12AspocOnR.x,...
      EfpiR.x-Es34AspocOffR.x,EfpiR.x-Es34AspocOnR.x},'comp')
  end
  ylabel(hca,'\Delta FPI x [mV/m]')
  set(hca,'YLim',[-4.9 .9])
  irf_legend(hca,{'12','12ao','34','34ao'},[0.95, 0.95])

  leg = {'12','12ao','34','34ao','dis'};
  hca = irf_panel('Ey'); set(hca,'ColorOrder',myCols)
  irf_plot(hca,{Es12AspocOffR.y,Es12AspocOnR.y,...
    Es34AspocOffR.y,Es34AspocOnR.y},'comp');
  if ~isempty(Efpi)
    hold(hca,'on'), irf_plot(hca,EfpiR.y,'m.'), hold(hca,'off'),
  end
  if ~isempty(Ehpca)
    hold(hca,'on'), irf_plot(hca,EhpcaR.y,'r.'), hold(hca,'off'), leg = [leg 'hpca'];
  end
  ylabel(hca,'Ey [mV/m]'), irf_legend(hca,leg,[0.95, 0.95])

  hca = irf_panel('Fpi_y'); set(hca,'ColorOrder',myCols)
  if ~isempty(Efpi)
    irf_plot(hca,{EfpiR.y-Es12AspocOffR.y,EfpiR.y-Es12AspocOnR.y,...
      EfpiR.y-Es34AspocOffR.y,EfpiR.y-Es34AspocOnR.y},'comp')
  end
  ylabel(hca,'\Delta FPI Y [mV/m]')
  set(hca,'YLim',[-1.9 1.9])
  irf_legend(hca,{'12','12ao','34','34ao'},[0.95, 0.95])

  %irf_plot_ylabels_align(h)
  irf_zoom(h,'x',Tint)
  mmsIdS = sprintf('MMS%d',mmsId);
  title(h(1),mmsIdS)

  set(gcf,'paperpositionmode','auto')
  irf_print_fig(['Offsets_' mmsIdS '_' irf_fname(Tint.start.epochUnix)],'png')
end
%% Scatter plots
if isempty(DeltaAspocOff)
  irf.log('warning','No data with ASPOC OFF')
  Es12AspocOffR = Es12AspocOnR; Es34AspocOffR = Es34AspocOnR;
  DeltaAspocOff = DeltaAspocOn;
end

Del = delta_off(DeltaAspocOff.data);

if ~any(idxMSH) || sum(idxMSH)<5 || all(all(isnan(Es12AspocOffR.data(idxMSH,:))))
  irf.log('warning','No MSH data')
  idxMSH = ~idxMSH;
end
TintTmp = Es12AspocOffR.time(idxMSH);
res = struct('tint',EpochUnix(median(TintTmp.epochUnix)),...
  'p12',[],'p34',[],'p1234',[],'delta',DeltaAspocOffRes,'ni', NifpiRes,...
  'psp',PSPRes);
h = irf_figure(94,4,'reset');
figure(94), clf
clf
nRows = 4;
subplot(nRows,2,1)
[~,res.p12.x] = plot_xy(Es12AspocOffR.x.data(idxMSH)-Del.x,EfpiR.x.data(idxMSH));
title(mmsIdS),ylabel('Ex FPI [mV/m]'), xlabel('SDP 12 [mV/m]')

subplot(nRows,2,2)
[~,res.p12.y] = plot_xy(Es12AspocOffR.y.data(idxMSH)-Del.y,EfpiR.y.data(idxMSH));
ylabel('Ey FPI [mV/m]'), xlabel('SDP 12 [mV/m]')
title([Tint.start.utc(1) ' - ' Tint.stop.utc(1)])

subplot(nRows,2,3)
[~,res.p34.x] = plot_xy(Es34AspocOffR.x.data(idxMSH),EfpiR.x.data(idxMSH));
ylabel('Ex FPI [mV/m]'), xlabel('SDP 34 [mV/m]')

subplot(nRows,2,4)
[~,res.p34.y] = plot_xy(Es34AspocOffR.y.data(idxMSH),EfpiR.y.data(idxMSH));
ylabel('Ey FPI [mV/m]'), xlabel('SDP 34 [mV/m]')

subplot(nRows,2,5)
[~,res.p12.xy,residuals] = plot_mvregress(...
  Es12AspocOffR.data(idxMSH,:)-repmat([Del.x Del.y],size(Es12AspocOffR.data(idxMSH,1)),1),EfpiR.data(idxMSH,:));
ylabel('E FPI [mV/m]'), xlabel('SDP 12 [mV/m]')
%Es12AspocOffResidual1 = irf.ts_vec_xy(Es12AspocOffR.time(idxMSH),residuals(:,1:2));

subplot(nRows,2,6)
[~,res.p34.xy,residuals] = plot_mvregress(...
  Es34AspocOffR.data(idxMSH,:),EfpiR.data(idxMSH,:));
ylabel('E FPI [mV/m]'), xlabel('SDP 34 [mV/m]')
%Es34AspocOffResidual1 = irf.ts_vec_xy(Es34AspocOffR.time(idxMSH),residuals(:,1:2));

if 0
  subplot(4,2,7) %#ok<UNRCH>
  [~,res.p1234,residuals] = plot_mvregress(...
    [Es12AspocOffR.data(idxMSH,:) Es34AspocOffR.data(idxMSH,:)],...
    [EfpiR.data(idxMSH,1:2) EfpiR.data(idxMSH,1:2)]);
  ylabel('E FPI [mV/m]'), xlabel('SDP [mV/m]')
  Es12AspocOffResidual2 = irf.ts_vec_xy(Es12AspocOffR.time(idxMSH),residuals(:,1:2));
  Es34AspocOffResidual2 = irf.ts_vec_xy(Es12AspocOffR.time(idxMSH),residuals(:,3:4));
end

if 1
  subplot(nRows,2,7)
  [~,res.p1234.x] = plot_xy([Es12AspocOffR.data(idxMSH,1)-Del.x; Es34AspocOffR.data(idxMSH,1)],...
    [EfpiR.x.data(idxMSH); EfpiR.x.data(idxMSH)]); cla
  [offs_x, slope] = comp_off_slope(Es12AspocOffR(idxMSH)-[Del.x Del.y], Es34AspocOffR(idxMSH), EfpiR(idxMSH));
  res.p1234.xy.slope = slope; res.p1234.xy.offs=offs_x;
  ylabel('E FPI [mV/m]'), xlabel('SDP [mV/m]')

  subplot(nRows,2,8)
  [~,res.p1234.y] = plot_xy([Es12AspocOffR.data(idxMSH,2)-Del.y; Es34AspocOffR.data(idxMSH,2)],...
    [EfpiR.y.data(idxMSH); EfpiR.y.data(idxMSH)]); cla
  if isempty(EhpcaR) || sum(~isnan(EhpcaR.data(idxMSH,1))) < 5, offs_x =[]; slope = [];
  else
    [offs_x, slope] = comp_off_slope(Es12AspocOffR(idxMSH)-[Del.x Del.y], Es34AspocOffR(idxMSH), EhpcaR(idxMSH));
  end
  res.p1234.xy.slope_hpca = slope; res.p1234.xy.offs_hpca=offs_x;
  ylabel('E HPCA H+ [mV/m]'), xlabel('SDP [mV/m]')
end

set(gcf,'paperpositionmode','auto')
irf_print_fig(['ScatterPlot_' mmsIdS '_' irf_fname(Tint.start.epochUnix)],'png')

%% Validation figure
% resample to 5 sec
%epoch5sec = fix(Tint.start.epochUnix/60)*60:5:ceil(Tint.stop.epochUnix/60)*60;
%Epoch5sec = EpochUnix(epoch5sec);
%EfpiR5 = Efpi.resample(Epoch5sec,'median');
if ~isempty(Es12AspocOff)
  Es12Corr = Es12AspocOff; Es34Corr = Es34AspocOff;
else
  Es12Corr = Es12AspocOn; Es34Corr = Es34AspocOn;
end
if 0
  Es12Corr.data = Es12Corr.data*res.p1234.slope; %#ok<UNRCH>
  Es12Corr.data(:,1) = Es12Corr.data(:,1) + res.p1234.offs(1);
  Es12Corr.data(:,2) = Es12Corr.data(:,2) + res.p1234.offs(2);
  Es34Corr.data = Es34Corr.data*res.p1234.slope;
  Es34Corr.data(:,1) = Es34Corr.data(:,1) + res.p1234.offs(3);
  Es34Corr.data(:,2) = Es34Corr.data(:,2) + res.p1234.offs(4);
  Es12AspocOffResidual = Es12AspocOffResidual2;
  Es34AspocOffResidual = Es34AspocOffResidual2;
else
  %alpha = res.p34.xy.slope;
  %OffDSL.x = - res.p34.xy.offs(1); OffDSL.y = - res.p34.xy.offs(2);
  alpha = res.p1234.xy.slope;
  OffDSL.x = res.p1234.xy.offs(1); OffDSL.y = 0;
  Es12Corr.data = Es12Corr.data*alpha;
  Es12Corr.data(:,1) = Es12Corr.data(:,1) -Del.x -OffDSL.x;
  Es12Corr.data(:,2) = Es12Corr.data(:,2) -Del.y -OffDSL.y;
  Es34Corr.data = Es34Corr.data*alpha;
  Es34Corr.data(:,1) = Es34Corr.data(:,1) -OffDSL.x;
  Es34Corr.data(:,2) = Es34Corr.data(:,2) -OffDSL.y;
  EfpiTmp = Efpi.resample(Es12Corr,'spline');
  Es12Residual = Es12Corr;
  Es12Residual.data = EfpiTmp.data(:,1:2) - Es12Corr.data;
  EfpiTmp = Efpi.resample(Es34Corr,'spline');
  Es34Residual = Es34Corr;
  Es34Residual.data = EfpiTmp.data(:,1:2) - Es34Corr.data;
  if ~isempty(Ehpca)
    EhpcaTmp = Ehpca.resample(Es12Corr,'spline');
    Es12ResidualHPCA = Es12Corr;
    Es12ResidualHPCA.data = EhpcaTmp.data(:,1:2) - Es12Corr.data;
    EhpcaTmp = Ehpca.resample(Es34Corr,'spline');
    Es34ResidualHPCA = Es34Corr;
    Es34ResidualHPCA.data = EhpcaTmp.data(:,1:2) - Es34Corr.data;
  end
end

h = irf_figure(95,8,'reset');
cmap = irf_colormap('space'); colormap(cmap)

hca = irf_panel('B');
irf_plot(hca,B);
if any(B.abs.data>100), set(hca,'YLim',[-99 99]), end
ylabel(hca,'B [nT]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

hca = irf_panel('ISpec');
if ~isempty(IonSpc)
  irf_spectrogram(hca,IonSpc,'log','donotfitcolorbarlabel');
end
set(hca,'YScale','log')
hold(hca,'on')
irf_plot(hca,Tifpi)
hold(hca,'off')
ylabel(hca,'E(i) [eV]'), set(hca,'YTick',[10 100 1000 10000])
irf_plot_axis_align(h);

hca = irf_panel('Vfpi');
irf_plot(hca,Vifpi);
ylabel(hca,'Vi [km/s]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

hca = irf_panel('Ni');
irf_plot(hca,Nifpi);
ylabel(hca,'Ni [cc]')
% ScPot
hca = axes('Position',get(hca,'Position'));
hl = irf_plot(hca,PSP); set(hl,'Color','b')
ylabel(hca,'PSP [V]')
set(hca, 'YAxisLocation','right');
set(hca,'Color','none'); % color of axis
set(hca,'XColor','b','YColor','b');
irf_zoom(hca,'x',Tint), xlabel(hca,'')
set(hca,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required

leg = {'FPI'}; co = [[1 0 1];[0 0 0];[0 0 1];[1 0 0]];
hca = irf_panel('Ex');
irf_plot(hca,Efpi.x,'m.')
hold(hca,'on')
if ~isempty(Ehpca),irf_plot(hca,EhpcaR.x,'r.'), leg = [leg 'hpca'];
  co = [[1 0 1];[1 0 0];[0 0 0];[0 0 1];];
end
leg = [leg, '12' , '34'];
irf_plot(hca,{Es12Corr.x,Es34Corr.x},'comp')
hold(hca,'off')
ylabel(hca,'Ex [mV/m]')
set(hca,'Ylim',9*[-1 1])
set(hca,'ColorOrder',co);
irf_legend(hca,leg,[0.95, 0.95])
irf_legend(hca,sprintf('\\alpha=%.2f, \\Delta=%.2f mV/m, DSL=%.2f mV/m',...
  alpha, Del.x, OffDSL.x),[0.05, 0.95])

%myCols = [[0 0 1];[1 0 0];];
hca = irf_panel('Ex_res'); %set(hca,'ColorOrder',myCols)
if ~isempty(Ehpca)
  irf_plot(hca,{Es12Residual.x,Es34Residual.x,...
    Es12ResidualHPCA.x,Es34ResidualHPCA.x},'comp')
else, irf_plot(hca,{Es12Residual.x,Es34Residual.x},'comp')
end
ylabel(hca,'Ex-res [mV/m]')
set(hca,'Ylim',1.99*[-1 1])
irf_legend(hca,{'12-dis','34-dis','12-hpca','34-hpca'},[0.95, 0.95])

leg = {'FPI'}; co = [[1 0 1];[0 0 0];[0 0 1];[1 0 0]];
hca = irf_panel('Ey');
irf_plot(hca,Efpi.y,'m.')
hold(hca,'on')
if ~isempty(Ehpca),irf_plot(hca,EhpcaR.y,'r.'), leg = [leg 'hpca'];
  co = [[1 0 1];[1 0 0];[0 0 0];[0 0 1];];
end
leg = [leg, '12' , '34'];
irf_plot(hca,{Es12Corr.y,Es34Corr.y},'comp')
hold(hca,'off')
ylabel(hca,'Ey [mV/m]')
set(hca,'Ylim',9*[-1 1])
set(hca,'ColorOrder',co);
irf_legend(hca,leg,[0.95, 0.95])
irf_legend(hca,sprintf('\\alpha=%.2f, \\Delta=%.2f mV/m, DSL=%.2f mV/m',...
  alpha, Del.y, OffDSL.y),[0.05, 0.95])

hca = irf_panel('Ey_res');  %set(hca,'ColorOrder',myCols)
if ~isempty(Ehpca)
  irf_plot(hca,{Es12Residual.y,Es34Residual.y,...
    Es12ResidualHPCA.y,Es34ResidualHPCA.y},'comp')
else, irf_plot(hca,{Es12Residual.y,Es34Residual.y},'comp')
end
ylabel(hca,'Ey-res [mV/m]')
set(hca,'Ylim',1.99*[-1 1])
irf_legend(hca,{'12-dis','34-dis','12-hpca','34-hpca'},[0.95, 0.95])

irf_zoom(h,'x',Tint)

title(h(1),mmsIdS)

set(gcf,'paperpositionmode','auto')
irf_print_fig(['Verification_' mmsIdS '_' irf_fname(Tint.start.epochUnix)],'png')
end



function ints = find_on(mask)

idxJump = find(diff(mask)~=0);

ints = []; iStop = [];
for i=1:length(idxJump)+1
  if isempty(iStop), iStart = 1; else, iStart = iStop + 1; end
  if i==length(idxJump)+1, iStop = length(mask); else, iStop = idxJump(i); end
  if ~mask(iStart), continue, end
  ints = [ ints; iStart iStop]; %#ok<AGROW>
end

end

function [h,fitParams] = plot_xy(x,y)

idx = isnan(x); x(idx) = []; y(idx) = [];
idx = isnan(y); x(idx) = []; y(idx) = [];

p=polyfit( x,y,1);
d = polyval(p,x)-y;
idx = abs(d)<3*std(d); idxOut = abs(d)>=3*std(d); % remove outlyers
p=polyfit( x(idx),y(idx),1);
cc=corrcoef(y(idx),x(idx));
slope = p(1); offs = p(2); corr_coef = cc(1,2);

h = plot(x(idx),y(idx),'.',x(idxOut),y(idxOut),'o');
ax=axis; set(gca,'YLim',10*[-1 1],'XLim',10*[-1 1])
ymax=ax(4);ymin=ax(3);dy=(ymax-ymin)/20;
ytext=ymax-dy;
xtext=ax(1)+(ax(2)-ax(1))/20;

hold on
xp=[min(x) max(x)];
plot(xp,polyval(p,xp),'k-');
axis(ax);
text(xtext,ytext,...
  ['slope=' num2str(slope,3) '  offs=' num2str(offs,2)]);ytext=ytext-dy;
text(xtext,ytext,['cc=' num2str(corr_coef,3)]);

fitParams = struct('slope',slope,'offs',offs,'corrCoef',corr_coef,'range',xp);
end

function [h,fitParams,residuals] = plot_mvregress(x,y)

residuals = nan(size(x));
idx = ~isnan(x(:,1)) & ~isnan(y(:,1));

x = x(idx,:); y = y(idx,:);

n=size(x,1);
d=size(x,2);
X = cell(n,1);

for i=1:n
  X{i} = [eye(d) x(i,1:d)'];
end
[beta,~,residualsTmp] = mvregress(X,y(:,1:d));
residuals(idx,:) = residualsTmp;

idxTmp = [];
if 1 % disregard points with > 3*std
  ia = abs(residualsTmp)<3*repmat(std(residualsTmp),n,1);
  idxTmp = ia(:,1);
  for i=2:d, idxTmp = idxTmp & ia(:,i); end
  [beta,~,residualsTmp2] = mvregress(X(idxTmp),y(idxTmp,1:d));
  residualsTmp(idxTmp,:) = residualsTmp2;
  residualsTmp(~idxTmp,:) = NaN;
  residuals(idx,:) = residualsTmp;
end

slope = beta(end); offs = beta(1:end-1); corr_coef = [];

%%
cols = 'krgb';
ax=axis; set(gca,'YLim',10*[-1 1],'XLim',10*[-1 1])
ax=axis;
ymax=ax(4);ymin=ax(3);dy=(ymax-ymin)/10;
ytext=ymax-dy; xtext=ax(1)+(ax(2)-ax(1))/20;

hold on

xp = zeros(d,2);
for i=1:d
  plot(x(:,i),y(:,i),[cols(i) '.'])
  if ~isempty(idxTmp), plot(x(~idxTmp,i),y(~idxTmp,i),[cols(i) 'o']), end
  if i==1
    text(xtext,ytext,['slope=' num2str(slope,3)]), ytext=ytext-dy;
  end
  text(xtext,ytext,['offs=' num2str(offs(i),2)]), ytext=ytext-dy;
  p = [slope offs(i)]; xp(i,:)=[min(x(:,i)) max(x(:,i))];
  plot(xp(i,:),polyval(p,xp(i,:)),[cols(i) '-']);
end
hold off
h = ax;
fitParams = struct('slope',slope,'offs',offs,'corrCoef',corr_coef,'range',xp);
end
