function mms_sdp_pl_delta_off(fileName)

d = dataobj(fileName);
mmsId = str2double(d.GlobalAttributes.Source_name{:}(4));

%%
bitmask = d.data.(sprintf('mms%d_edp_dce_bitmask',mmsId)).data;
mskAsponOn = bitand(bitmask(:,2),64) > 0;
mskAsponOn = mskAsponOn | (bitand(bitmask(:,3),64) > 0);

epochSpin = d.data.(sprintf('mms%d_edp_dce_spinfit_epoch',mmsId)).data;
epochFull = d.data.(sprintf('mms%d_edp_dce_epoch',mmsId)).data;
Tint = EpochTT(epochFull([1 end]));

p12 =  d.data.(sprintf('mms%d_edp_dce_spinfit_e12',mmsId)).data;
p12(p12<-1e20) = NaN;
p34 =  d.data.(sprintf('mms%d_edp_dce_spinfit_e34',mmsId)).data;
p34(p34<-1e20) = NaN;

%% Resample bitmask
spinSize = int64(20000000000);
mskAsponOnSpin = zeros(size(epochSpin));
ints = find_on(mskAsponOn);
for i=1:size(ints,1)
  mskAsponOnSpin((epochSpin> epochFull(ints(1,1))-spinSize/2) & ... 
    (epochSpin<= epochFull(ints(1,2))+spinSize/2)) = 1;
end
mskAsponOnSpin=logical(mskAsponOnSpin);

%%
EpochS = EpochTT(epochSpin);
Es12 = irf.ts_vec_xy(EpochS,p12(:,3:4));
Es34 = irf.ts_vec_xy(EpochS,p34(:,3:4));
Es12AspocOff =  Es12(~mskAsponOnSpin); Es12AspocOn =  Es12(mskAsponOnSpin);
Es34AspocOff =  Es34(~mskAsponOnSpin); Es34AspocOn =  Es34(mskAsponOnSpin);
DeltaAspocOff = Es12AspocOff - Es34AspocOff;
DeltaAspocOn = Es12AspocOn - Es34AspocOn;

%% FPI
B = mms.get_data('dfg_ql_srvy',Tint,mmsId);
Vfpi = mms.get_data('Vi_gse_fpi_ql',Tint,mmsId);
if isempty(Vfpi)
  Nifpi = []; Efpi = [];
else
  Nifpi = mms.get_data('Ni_fpi_ql',Tint,mmsId);
  Efpi = irf_e_vxb(Vfpi,B.resample(Vfpi));
end

%% resample to 1 min
epoch1min = fix(Tint.start.epochUnix/60):ceil(Tint.stop.epochUnix/60);
epoch1min = epoch1min*60;
Epoch1min = EpochUnix(epoch1min);


if ~isempty(Es12AspocOff), Es12AspocOffR = Es12AspocOff.resample(Epoch1min,'median');
else Es12AspocOffR = Es12AspocOff;
end
if ~isempty(Es12AspocOn), Es12AspocOnR = Es12AspocOn.resample(Epoch1min,'median');
else Es12AspocOnR = Es12AspocOn;
end
if ~isempty(Es34AspocOff), Es34AspocOffR = Es34AspocOff.resample(Epoch1min,'median');
else Es34AspocOffR = Es34AspocOff;
end
if ~isempty(Es34AspocOn), Es34AspocOnR = Es34AspocOn.resample(Epoch1min,'median');
else Es34AspocOnR = Es34AspocOn;
end
if ~isempty(Efpi)
  EfpiR = Efpi.resample(Epoch1min,'median');
end
if isempty(DeltaAspocOff), DeltaAspocOffR = DeltaAspocOff;
else DeltaAspocOffR = DeltaAspocOff.resample(Epoch1min,'median');
end
if isempty(DeltaAspocOn), DeltaAspocOnR = DeltaAspocOn;
else DeltaAspocOnR = DeltaAspocOn.resample(Epoch1min,'median');
end

%%
myCols = [[0 0 0];[.3 .3 .3];[0 0 1];[.2 .2 .8]];
h = irf_figure(93,6);
hca = irf_panel('B');
irf_plot(hca,B);
if any(B.abs.data>100), set(hca,'YLim',[-99 99]), end
ylabel(hca,'B [nT]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

hca = irf_panel('Vfpi');
irf_plot(hca,Vfpi);
ylabel(hca,'Vi [km/s]'), irf_legend(hca,{'X','Y','Z'},[0.95, 0.95])

hca = irf_panel('Ni');
irf_plot(hca,Nifpi);
ylabel(hca,'Ni [cc]')

hca = irf_panel('Delta_off'); set(hca,'ColorOrder',myCols)
irf_plot(hca,{DeltaAspocOffR.x,DeltaAspocOnR.x,DeltaAspocOffR.y,DeltaAspocOnR.y},'comp')
ylabel(hca,'\Delta Off [mV/m]')
set(hca,'YLim',[-1.4 1.4])
irf_legend(hca,{'X','Y'},[0.95, 0.95])

hca = irf_panel('Fpi_x'); set(hca,'ColorOrder',myCols)
if ~isempty(Efpi)
  irf_plot(hca,{EfpiR.x-Es12AspocOffR.x,EfpiR.x-Es12AspocOnR.x,...
    EfpiR.x-Es34AspocOffR.x,EfpiR.x-Es34AspocOnR.x},'comp')
end
ylabel(hca,'\Delta FPI x [mV/m]')
set(hca,'YLim',[-4.9 .9])
irf_legend(hca,{'12','34'},[0.95, 0.95])

hca = irf_panel('Fpi_y'); set(hca,'ColorOrder',myCols)
if ~isempty(Efpi)
  irf_plot(hca,{EfpiR.y-Es12AspocOffR.y,EfpiR.y-Es12AspocOnR.y,...
    EfpiR.y-Es34AspocOffR.y,EfpiR.y-Es34AspocOnR.y},'comp')
end
ylabel(hca,'\Delta FPI Y [mV/m]')
set(hca,'YLim',[-1.9 1.9])
irf_legend(hca,{'12','34'},[0.95, 0.95])

%irf_plot_ylabels_align(h)
irf_zoom(h,'x',Tint)
mmsIdS = sprintf('MMS%d',mmsId);
title(h(1),mmsIdS)

set(gcf,'paperpositionmode','auto')
print('-dpng',['DeltaOff_' mmsIdS '_' irf_fname(Tint.start.epochUnix)])

%%

NifpiR = Nifpi.resample(Epoch1min,'median');
idxHighN = NifpiR.data>5;

end

function ints = find_on(mask)

%%

idxJump = find(diff(mask)~=0);

ints = []; iStop = [];
for i=1:length(idxJump)+1
  if isempty(iStop), iStart = 1; else iStart = iStop + 1; end
  if i==length(idxJump)+1, iStop = length(mask); else iStop = idxJump(i); end
  if ~mask(iStart), continue, end
  ints = [ ints; iStart iStop]; %#ok<AGROW>
end

end
