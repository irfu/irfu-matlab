function res = mms_sdp_comp_msh_dsl_off(Tint)
%MMS_SDP_COMP_MSH_DSL_OFF  DSL ofssets from FPI comparison
%
% res = MMS_SDP_COMP_MSH_DSL_OFF(Tint)

res = struct('c1',[],'c2',[],'c3',[],'c4',[],'tint',Tint);

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

%%
epoch1min = fix(Tint.start.epochUnix/60)*60:20:ceil(Tint.stop.epochUnix/60)*60;
%epoch1min = epoch1min*60;
Epoch20s = EpochUnix(epoch1min);

idxMSH = [];
E34 = struct('c1',[],'c2',[],'c3',[],'c4',[]);
EFPI = E34; flagOldFile = false;

maneuvers = mms_maneuvers(Tint); % Locate any maneuvers

% Check if Tint covers midnight and the stop time is more than 10 minutes
% after midnight and within the last 40 days. If this is the case then
% loading recent L2 or L2Pre data products may result in data only from the
% start of Tint to midnight if the next day(-s) L2 or L2Pre products has
% not yet been produced.
if( ~strcmp(Tint.start.utc('dd'), Tint.stop.utc('dd')) && ...
    EpochTT(Tint.stop.utc('yyyy-mm-ddT00:10:00.000000000Z')) < Tint.stop && ...
    Tint.stop.ttns > irf_time(now,'datenum>ttns') - int64(40*86400*10^9) )
  midnightROI = true;
else
  midnightROI = false;
end

for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  fPre = sprintf('mms%d_edp_fast_l2a_dce2d',mmsId);

  if flagOldFile
    Bmask = mms.db_get_ts(fPre,...
      sprintf('mms%d_edp_dce_bitmask',mmsId), Tint);
  else
    Bmask = mms.db_get_ts(fPre,...
      sprintf('mms%d_edp_bitmask_fast_l2a',mmsId), Tint);
    if mmsId==1 && isempty(Bmask)
      Bmask = mms.db_get_ts(fPre,...
        sprintf('mms%d_edp_dce_bitmask',mmsId), Tint);
      if ~isempty(Bmask)
        irf.log('warning','Using old L2a file')
        flagOldFile = true;
      end
    end
  end

  if ~isempty(Bmask)
    if flagOldFile
      P34 = mms.db_get_ts(fPre,...
        sprintf('mms%d_edp_dce_spinfit_e34',mmsId), Tint);
      p34 = P34.data(:,3:4);
    else
      P34 = mms.db_get_ts(fPre,...
        sprintf('mms%d_edp_espin_p34_fast_l2a',mmsId), Tint);
      p34 = P34.data;
    end
    epochSpin = P34.time.ttns;
    bitmask = Bmask.data;
    mskAsponOn = bitand(bitmask(:,2), MMS_CONST.Bitmask.ASPOC_RUNNING) > 0;
    mskAsponOn = mskAsponOn | (bitand(bitmask(:,3), MMS_CONST.Bitmask.ASPOC_RUNNING) > 0);
    epochFull = Bmask.time.ttns;
    % XXX hack
    if ~any(mskAsponOn)
      TAspOn = Tint.stop+(-3*3600); mskAsponOn = epochFull>TAspOn.ttns;
    end

    % Resample bitmask
    spinSize = MMS_CONST.Limit.SPINFIT_INTERV;
    mskAsponOnSpin = false(size(epochSpin));
    ints = find_on(mskAsponOn);
    for i=1:size(ints,1)
      mskAsponOnSpin((epochSpin> epochFull(ints(1,1))-spinSize/2) & ...
        (epochSpin<= epochFull(ints(1,2))+spinSize/2)) = true;
    end
    EpochS = EpochTT(epochSpin);
    Es34 = irf.ts_vec_xy(EpochS,p34(:,1:2));
    Es34AspocOff =  Es34(~mskAsponOnSpin);
    if isempty(Es34AspocOff) || all(all(isnan(Es34AspocOff.data)))
      irf.log('warning','No data w/o ASPOC, using all data')
      Es34AspocOff = Es34;
    end
    if ~isempty(Es34AspocOff)
      Es34AspocOffR = Es34AspocOff.resample(Epoch20s,'median');
    else, Es34AspocOffR = [];
    end
    E34.(mmsIdS) = Es34AspocOffR;
  end

  for iManuev = 1:length(maneuvers.(sprintf('mms%i',mmsId)))
    % Maneuvers was found...
    irf.log('warning', ['Maneuver: ', ...
      maneuvers.(sprintf('mms%i',mmsId)){iManuev}.start.toUtc, '/', ...
      maneuvers.(sprintf('mms%i',mmsId)){iManuev}.stop.toUtc]);
    % Set EDP data as NaN for the duration of the maneuvers.
    E34.(mmsIdS).data(bitand( E34.(mmsIdS).time.ttns >= maneuvers.(sprintf('mms%i',mmsId)){iManuev}.start.ttns, ...
      E34.(mmsIdS).time.ttns <= maneuvers.(sprintf('mms%i',mmsId)){iManuev}.stop.ttns), :) = NaN;
  end

  B = mms.get_data('B_dmpa_srvy',Tint,mmsId);
  if isempty(B) || (midnightROI && B.time(end) < EpochTT(Tint.stop.utc('yyyy-mm-ddT00:01:00.000000000Z')))
    % If higher level data was not found, or ROI covers midnight +10min and
    % resulting data did not at least cover midnight +1min it is likely
    % that the higher level has not yet been produced. Try a lower level
    % product.
    B = mms.get_data('B_dmpa_dfg_srvy_l2pre',Tint,mmsId);
    if isempty(B) || (midnightROI && B.time(end) < EpochTT(Tint.stop.utc('yyyy-mm-ddT00:01:00.000000000Z')))
      B = mms.get_data('B_dmpa_dfg_srvy_ql',Tint,mmsId);
      if isempty(B)
        irf.log('warning','Failed to load any B-field data');
        continue
      else
        irf.log('warning','Using QL DFG data');
      end
    else
      irf.log('warning','Using L2pre DFG data');
    end
  end

  % Here we try different sources of FPI data
  Vifpi = mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,mmsId);
  if isempty(Vifpi) || (midnightROI && Vifpi.time(end) < EpochTT(Tint.stop.utc('yyyy-mm-ddT00:01:00.000000000Z')))
    % No L2, try L1b
    %Vifpi = mms.get_data('Vi_gse_fpi_fast_l1b',Tint,mmsId);
    %if ~isempty(Vifpi)
    %  irf.log('warning','Using L1b FPI data');
    %  if isempty(idxMSH), Nifpi = mms.get_data('Ni_fpi_fast_l1b',Tint,mmsId); end
    %else
    % No L2 or L1b, try QL
    Vifpi = mms.get_data('Vi_dbcs_fpi_ql',Tint,mmsId);
    if ~isempty(Vifpi)
      irf.log('warning','Using QL FPI data');
      if isempty(idxMSH), Nifpi = mms.get_data('Ni_fpi_ql',Tint,mmsId); end
    else
      % No L2, L1b or QL. Last resort, try SITL
      %Vifpi = mms.get_data('Vi_dbcs_fpi_sitl',Tint,mmsId);
      %if ~isempty(Vifpi)
      %irf.log('warning','Using SITL FPI data');
      %if isempty(idxMSH), Nifpi = mms.get_data('Ni_fpi_sitl',Tint,mmsId); end
      %else
      irf.log('warning', 'Did not find any FPI data');
      continue
    end
    %end
    %end
  else
    irf.log('notice','Using L2 FPI data');
    if isempty(idxMSH), Nifpi = mms.get_data('Ni_fpi_fast_l2',Tint,mmsId); end
  end
  Efpi = irf_e_vxb(Vifpi,B.resample(Vifpi));
  EfpiR = Efpi.resample(Epoch20s,'median');

  EFPI.(mmsIdS) = EfpiR;

  if isempty(idxMSH)
    NifpiR = Nifpi.resample(Epoch20s,'median');
    VifpiR = Vifpi.resample(Epoch20s,'median');
    idxMSH = NifpiR.data>5 & VifpiR.x.data>-200;
    if sum(idxMSH) < 180 % We require min 100 spins of MSH data
      irf.log('warning','Using ALL data (not only MSH)')
      idxMSH(:)=true;
    end
  end
end
%%
if isempty(idxMSH), return, end

ErefFpiX = []; ErefFpiY = []; nData = [];
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(EFPI.(mmsIdS)), continue, end
  if isempty(nData), nData = size(EFPI.(mmsIdS).data,1);
  elseif nData ~= size(EFPI.(mmsIdS).data,1)
    msgS = ['different number of point in FPI on ' mmsIdS];
    irf.log('critical',msgS), error(msgS)
  end
  ErefFpiX = [ErefFpiX EFPI.(mmsIdS).x.data]; %#ok<AGROW>
  ErefFpiY = [ErefFpiY EFPI.(mmsIdS).y.data]; %#ok<AGROW>
end

ErefFpiX = median(ErefFpiX,2,'omitnan'); ErefFpiY = median(ErefFpiY,2,'omitnan');
ErefFpi = irf.ts_vec_xy(Epoch20s,[ErefFpiX ErefFpiY]);

%%
if flagOldFile, ALPHA = 1.25/1.1; else, ALPHA = 1; end
DE = struct('c1',[],'c2',[],'c3',[],'c4',[]);
Off = DE;
for mmsId = 1:4
  mmsIdS = sprintf('c%d',mmsId);
  if isempty(E34.(mmsIdS)), continue, end
  DE.(mmsIdS) = ALPHA*E34.(mmsIdS) - ErefFpi;
  deTmp = DE.(mmsIdS).data; deTmp = deTmp(idxMSH,:);
  offTmp = zeros(1,2);
  for iComp = 1:2
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
  hold(hca,'on'), irf_plot(hca,ErefFpi.x,'c.'), hold(hca,'off')
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
  hold(hca,'on'), irf_plot(hca,ErefFpi.y,'c.'), hold(hca,'off')
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
    if isempty(Off.(mmsIdS)), tt = NaN; else, tt = Off.(mmsIdS)(1); end
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
    if isempty(Off.(mmsIdS)), tt = NaN; else, tt = Off.(mmsIdS)(2); end
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
  hold(hca,'on'), irf_plot(hca,ttt,'c.'), irf_plot(hca,ttt(idxMSH),'m.'), hold(hca,'off')
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
  hold(hca,'on'), irf_plot(hca,ttt,'c.'), irf_plot(hca,ttt(idxMSH),'m.'), hold(hca,'off')
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
  if isempty(iStop), iStart = 1; else, iStart = iStop + 1; end
  if i==length(idxJump)+1, iStop = length(mask); else, iStop = idxJump(i); end
  if ~mask(iStart), continue, end
  ints = [ ints; iStart iStop]; %#ok<AGROW>
end

end