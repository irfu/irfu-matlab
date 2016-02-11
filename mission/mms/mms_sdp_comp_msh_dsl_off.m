function res = mms_sdp_comp_msh_dsl_off(Tint)
%MMS_SDP_COMP_MSH_DSL_OFF  DSL ofssets from FPI comparison
%
% res = MMS_SDP_COMP_MSH_DSL_OFF(Tint,mmsId)

res = struct('c1',[],'c2',[],'c3',[],'c4',[],'tint',Tint);

%%
epoch1min = fix(Tint.start.epochUnix/60)*60:20:ceil(Tint.stop.epochUnix/60)*60;
%epoch1min = epoch1min*60;
Epoch20s = EpochUnix(epoch1min);

idxMSH = [];
E34 = struct('c1',[],'c2',[],'c3',[],'c4',[]);
EFPI = E34;

for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  fPre = sprintf('mms%d_edp_fast_l2a_dce2d',mmsId);
  Bmask = mms.db_get_ts(fPre,...
    sprintf('mms%d_edp_dce_bitmask',mmsId), Tint);
  if ~isempty(Bmask)
    P34 = mms.db_get_ts(fPre,...
      sprintf('mms%d_edp_dce_spinfit_e34',mmsId), Tint);
    p34 = P34.data;
    epochSpin = P34.time.ttns;
    bitmask = Bmask.data;
    mskAsponOn = bitand(bitmask(:,2),64) > 0;
    mskAsponOn = mskAsponOn | (bitand(bitmask(:,3),64) > 0);
    epochFull = Bmask.time.ttns;
    % XXX hack
    if ~any(mskAsponOn)
      TAspOn = Tint.stop+(-3*3600); mskAsponOn = epochFull>TAspOn.ttns;
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
    EpochS = EpochTT(epochSpin);
    Es34 = irf.ts_vec_xy(EpochS,p34(:,3:4));
    Es34AspocOff =  Es34(~mskAsponOnSpin);
    
    if ~isempty(Es34AspocOff),
      Es34AspocOffR = Es34AspocOff.resample(Epoch20s,'median');
    else Es34AspocOffR = [];
    end
    E34.(mmsIdS) = Es34AspocOffR;
  end
  
  B = mms.get_data('B_dmpa_srvy',Tint,mmsId);
  if isempty(B), B = mms.get_data('dfg_ql_srvy',Tint,mmsId); end
  if isempty(B), continue, end
  
  % XXX here we could try different sources of FPI data
  Vifpi = mms.get_data('Vi_gse_fpi_ql',Tint,mmsId);
  if isempty(Vifpi), continue, end
  Efpi = irf_e_vxb(Vifpi,B.resample(Vifpi));
  EfpiR = Efpi.resample(Epoch20s,'median');
  
  EFPI.(mmsIdS) = EfpiR;
  
  if isempty(idxMSH)
    Nifpi = mms.get_data('Ni_fpi_ql',Tint,mmsId);
    NifpiR = Nifpi.resample(Epoch20s,'median');
    VifpiR = Vifpi.resample(Epoch20s,'median');
    idxMSH = NifpiR.data>5 & VifpiR.x.data>-200;  
  end
end
%%
if isempty(idxMSH), return, end

ErefFpiX = []; ErefFpiY = [];
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(EFPI.(mmsIdS)), continue, end
  ErefFpiX = [ErefFpiX EFPI.(mmsIdS).x.data]; %#ok<AGROW>
  ErefFpiY = [ErefFpiY EFPI.(mmsIdS).y.data]; %#ok<AGROW>
end

ErefFpiX = median(ErefFpiX,2); ErefFpiY = median(ErefFpiY,2);
ErefFpi = irf.ts_vec_xy(Epoch20s,[ErefFpiX ErefFpiY]);

%%
ALPHA = 1.25/1.1;
DE = struct('c1',[],'c2',[],'c3',[],'c4',[]);
Off = DE;
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), continue, end
  DE.(mmsIdS) = ALPHA*E34.(mmsIdS) - ErefFpi;
  deTmp = DE.(mmsIdS).data; deTmp = deTmp(idxMSH,:); 
  offTmp = zeros(1,2);
  for iComp = 1:2,
    offTmp(iComp) = median(deTmp(~isnan(deTmp(:,iComp)),iComp));
  end
  Off.(mmsIdS) = offTmp;
end
res = Off; res.tint = Tint;

%% Validation figure
mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
h = irf_figure(195,6,'reset');

if 1 % Panel Ex
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {E34.(mmsIdS).x}]; %#ok<AGROW>
end

hca = irf_panel('Ex'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
hold(hca,'on'), irf_plot(hca,ErefFpi.x,'.'), hold(hca,'off')
ylabel(hca,'Ex [mV/m]')
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98, 0.1],'color','cluster');
end

if 1 % Panel Ey
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {E34.(mmsIdS).y}]; %#ok<AGROW>
end

hca = irf_panel('Ey'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
hold(hca,'on'), irf_plot(hca,ErefFpi.y,'.'), hold(hca,'off')
ylabel(hca,'Ey [mV/m]')
end

if 1 % Panel dEx
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(DE.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {DE.(mmsIdS).x}]; %#ok<AGROW>
end

hca = irf_panel('dEx'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
ylabel(hca,'dEx [mV/m]')

plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(Off.(mmsIdS)), tt = NaN; else tt = Off.(mmsIdS)(1); end
  plData= [plData {sprintf('dEx%d=%.2f',mmsId,tt)}]; %#ok<AGROW>
end
irf_legend(hca,plData,[0.1, 0.05],'color','cluster');
end

if 1 % Panel dEy
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(DE.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {DE.(mmsIdS).y}]; %#ok<AGROW>
end

hca = irf_panel('dEy'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
ylabel(hca,'dEy [mV/m]')

plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(Off.(mmsIdS)), tt = NaN; else tt = Off.(mmsIdS)(2); end
  plData= [plData {sprintf('dEy%d=%.2f',mmsId,tt)}]; %#ok<AGROW>
end
irf_legend(hca,plData,[0.1, 0.05],'color','cluster');
end

if 1 % Panel Ex-corr
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {E34.(mmsIdS).x*ALPHA-Off.(mmsIdS)(1)}]; %#ok<AGROW>
end

hca = irf_panel('Ex-corr'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
ttt = ErefFpi.x;
hold(hca,'on'), irf_plot(hca,ttt,'.'), irf_plot(hca,ttt(idxMSH),'.'), hold(hca,'off')
ylabel(hca,'Ex [mV/m]')
end

if 1 % Panel Ey-corr
plData = {};
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), plData= [plData {[]}]; continue, end %#ok<AGROW>
  plData= [plData {E34.(mmsIdS).y*ALPHA-Off.(mmsIdS)(2)}]; %#ok<AGROW>
end

hca = irf_panel('Ey-corr'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,plData,'comp')
ttt = ErefFpi.y;
hold(hca,'on'), irf_plot(hca,ttt,'.'), irf_plot(hca,ttt(idxMSH),'.'), hold(hca,'off')
ylabel(hca,'Ey [mV/m]')
end

irf_zoom(h,'x',Tint)

title(h(1),'MSH DSL Offsets')

set(gcf,'paperpositionmode','auto')
irf_print_fig(['SH_DSL_OFF_' irf_fname(Tint.start.epochUnix)],'png')

end

function ints = find_on(mask)

idxJump = find(diff(mask)~=0);

ints = []; iStop = [];
for i=1:length(idxJump)+1
  if isempty(iStop), iStart = 1; else iStart = iStop + 1; end
  if i==length(idxJump)+1, iStop = length(mask); else iStop = idxJump(i); end
  if ~mask(iStart), continue, end
  ints = [ ints; iStart iStop]; %#ok<AGROW>
end

end