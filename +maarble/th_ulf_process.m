function th_ulf_process(TT,thId,freqRange)
%TH_ULF_PROCESS  process THEMIS ULF data
%
%  th_ulf_process(TT,thId,freqRange)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% This software was developed as part of the MAARBLE (Monitoring,
% Analyzing and Assessing Radiation Belt Energization and Loss)
% collaborative research project which has received funding from the
% European Community's Seventh Framework Programme (FP7-SPACE-2011-1)
% under grant agreement n. 284520.



%% defaults
if nargin < 1
  %TT = [iso2epoch('2008-04-20T17:40:00Z') iso2epoch('2008-04-20T23:30:00Z')]; thId = 'e'; freqRange = 'pc35';
  TT = [iso2epoch('2008-04-20T22:20:00Z') iso2epoch('2008-04-20T23:30:00Z')]; thId = 'e'; freqRange = 'pc12';
elseif nargin < 3
  freqRange = 'all';
end

%% setup
dataDir = '/data/themis';
plotFlag = 1;
exportFlag = 1;

wantPC12 = 0;
wantPC35 = 0;
if ischar(freqRange)
  switch lower(freqRange)
    case 'all'
      wantPC12 = 1;
      wantPC35 = 1;
    case 'pc12'
      wantPC12 = 1;
    case 'pc35'
      wantPC35 = 1;
    otherwise
      error('Invalid value for freqRange')
  end
else
  error('frequency range must be one of: pc12, pc35, all')
end

tmpR = load(sprintf('%s%smRth.mat',dataDir,filesep), '-mat', ['Rth' thId]);
gseR = tmpR.(['Rth' thId]);
gseV=gseR(1:end-1,:); dtR = diff(gseR(:,1));
for iComp=2:4, gseV(:,iComp)=diff(gseR(:,iComp))./dtR; end
clear tmpR dtR

if ~isa(TT,'irf.TimeTable'), TT=irf.TimeTable(TT); end

for ievent=1:numel(TT)
  tint=TT.TimeInterval(ievent,:);
  irf.log('warning',sprintf('processing %s\n',irf_disp_iso_range(tint,1)))


  %% Load data
  axxx = ls('/data/caalocal'); axxx = ls('/data/themis');
  % Round time interval to minutes
  tint = [floor(tint(1)/60) ceil(tint(2)/60)]*60;
  % Extend time interval by these ranges to avoid edge effects
  DT_PC5 = 80*60; DT_PC2 = 120;

  bs = th_read_l2(['th' thId '_fgs_dsl'],tint+DT_PC5*[-1 1]);
  % Clean backward time jumps, example THD
  bs = clear_backward_jump(bs,'BS');
  [bs,~,ttGap] = th_clean_eb(bs);
  if isempty(bs)
    irf.log('warning','skipping, no BS data'),continue,
  end
  if wantPC35
    es = th_read_l2(['th' thId '_efs_dot0_dsl'],tint+DT_PC5*[-1 1]);
    es = clear_backward_jump(es,'ES');
    es=th_clean_eb(es,'NaN');
  end
  if wantPC12
    bl = th_read_l2(['th' thId '_fgl_dsl'],tint+DT_PC2*[-1 1]);
    bl = clear_backward_jump(bl,'BL');
    ef = th_read_l2(['th' thId '_eff_dot0_dsl'],tint+DT_PC2*[-1 1]);
    if isempty(ef)
      irf.log('warning','no EF data')
    else
      % Remove backward time jumps in EF
      % Example THE 2007-08-11T09:40:02.489975Z (3 points)
      ef = clear_backward_jump(ef,'EF');
    end
  end


  R = th_gse2dsl(irf_tlim(gseR,tint+DT_PC5*[-1 1]),thId);
  V = th_gse2dsl(irf_tlim(gseV,tint+DT_PC5*[-1 1]),thId);

  %% Calculate and plot
  bf = irf_filt(bs,0,1/600,1/5,5);
  t_1min = ((tint(1)-DT_PC5):60:(tint(end)+DT_PC5))';
  B0_1MIN = irf_resamp(bf,t_1min); %clear bf
  nGaps = size(ttGap.TimeInterval,1);
  if nGaps>0
    for iGap=1:nGaps
      tintGap = ttGap.TimeInterval(iGap,:);
      if diff(tintGap)>300
        irf.log('warning',['Long ' ttGap.Comment{iGap} ' in BS at '...
          irf_disp_iso_range(tintGap,1) ' > 5 min'])
        B0_1MIN = irf_tlim(B0_1MIN,tintGap+[-30 30],1);
      end
    end
  end
  if isempty(B0_1MIN), irf.log('warning','no BS data, skipping interval'), continue, end
  facMatrix = irf_convert_fac([],B0_1MIN,R);

  if wantPC35
    t_1SEC = ((tint(1)+2-DT_PC5):1:(tint(end)+DT_PC5))';
    B_1SEC = irf_resamp(bs,t_1SEC);

    % Set Gaps to NaN. for safety we set both E and B
    if nGaps>0
      for iGap=1:nGaps
        tintGap = ttGap.TimeInterval(iGap,:);
        B_1SEC(B_1SEC(:,1)>=tintGap(1)-0.5 & B_1SEC(:,1)<=tintGap(2)+0.5,2:end) = NaN;
        if ~isempty(es)
          es(es(:,1)>=tintGap(1)-0.5 & es(:,1)<=tintGap(2)+0.5,2:end) = NaN;
        end
      end
    end
    ii = find(diff(bs(:,1))>3*median(diff(bs(:,1))));
    if ~isempty(ii)
      for iGap=ii'
        irf.log('warning',['Long gap in BS ' irf_disp_iso_range(bs(iGap+[0 1],1)')])
        B_1SEC(t_1SEC(:,1)>bs(iGap,1) & t_1SEC(:,1)<bs(iGap+1,1),2:4) = NaN;
      end
    end
    B_1SEC(t_1SEC(:,1)<bs(1,1),2:4) = NaN; B_1SEC(t_1SEC(:,1)>bs(end,1),2:4) = NaN;

    if isempty(es), iE3D_1SEC = [];
    else
      E3D_1SEC = irf_resamp(es,t_1SEC);
      % Construct the inertial frame
      evxb = irf_tappl(irf_cross(B_1SEC,irf_resamp(V,t_1SEC)),'*1e-3*(-1)');
      iE3D_1SEC = E3D_1SEC;
      iE3D_1SEC(:,2:4) = iE3D_1SEC(:,2:4) - evxb(:,2:4);
      % Remove all E-fields > 10 mV/m in inertial frame
      for iComp=2:4, iE3D_1SEC(abs(iE3D_1SEC(:,iComp))>10, 2:4) = NaN; end
      ii = find(diff(es(:,1))>3*median(diff(es(:,1))));
      if ~isempty(ii)
        for iGap=ii'
          irf.log('warning',['Long gap in ES ' irf_disp_iso_range(es(iGap+[0 1],1)')])
          iE3D_1SEC(t_1SEC(:,1)>es(iGap,1) & t_1SEC(:,1)<es(iGap+1,1),2:4) = NaN;
        end
      end
      iE3D_1SEC(t_1SEC(:,1)<es(1,1),2:4) = NaN; iE3D_1SEC(t_1SEC(:,1)>es(end,1),2:4) = NaN;
    end

    ebsp = ...
      irf_ebsp(iE3D_1SEC,B_1SEC,[],B0_1MIN,[],'pc35',...
      'fac','polarization','noresamp','fullB=dB','facMatrix',facMatrix);
    ebsp.r = gseR; % add position for plotting
    tlim_ebsp(tint);
    if plotFlag
      close all
      h = irf_pl_ebsp(ebsp);
      irf_zoom(h,'x',tint)
      title(h(1),['THEMIS ' upper(thId) ', ' irf_disp_iso_range(tint,1)])
      orient tall
      print('-dpng',['MAARBLE_TH' upper(thId) '_ULF_PC35_' irf_fname(tint,5)])
    end
    if exportFlag
      maarble.export(ebsp,tint,['th' thId],'pc35')
    end
  end
  if wantPC12 && isempty(bl)
    irf.log('warning','skipping PC12, no BL data')
  elseif wantPC12 && ~isempty(bl)
    baseFreq = 16;
    fSampB = 1/median(diff(bl(:,1)));
    if isempty(ef)
      fSampE = NaN;
      t_BASE = (fix(bl(1,1)):1/baseFreq:ceil(bl(end,1)))';
      iE3D_BASE = [];
    else
      fSampE = 1/median(diff(ef(:,1)));
      t_BASE = (fix(min(bl(1,1),ef(1,1))):1/baseFreq:ceil(max(bl(end,1),ef(end,1))))';
    end

    B_BASE = irf_resamp(bl,t_BASE);
    T_TMP = t_BASE(t_BASE(:,1)>=bl(1,1) & t_BASE(:,1)<=bl(end,1)); ttB = T_TMP([1 end],1)';
    ii = find(diff(bl(:,1))>2/fSampB);
    if ~isempty(ii)
      for iGap=ii'
        irf.log('warning',['Long gap in BL ' irf_disp_iso_range(bl(iGap+[0 1],1)')])
        B_BASE(t_BASE(:,1)>bl(iGap,1) & t_BASE(:,1)<bl(iGap+1,1),2:4) = NaN;
        if bl(iGap+1,1)-bl(iGap,1)>DT_PC2
          irf.log('warning','splitting BL')
          T_TMP = t_BASE(t_BASE(:,1)>=bl(iGap,1) & t_BASE(:,1)<=bl(iGap+1,1));
          tint_TMP = [T_TMP(end) ttB(end,2)];
          ttB(end,2) = T_TMP(1); ttB = [ttB; tint_TMP];
        end
      end
    end
    B_BASE(t_BASE(:,1)<bl(1,1),2:4) = NaN; B_BASE(t_BASE(:,1)>bl(end,1),2:4) = NaN;

    ttE = [];
    if ~isempty(ef)
      E3D_BASE = irf_resamp(ef,t_BASE);
      T_TMP = t_BASE(t_BASE(:,1)>=ef(1,1) & t_BASE(:,1)<=ef(end,1)); ttE = T_TMP([1 end],1)';
      ii = find(diff(ef(:,1))>2/fSampE);
      if ~isempty(ii)
        for iGap=ii'
          irf.log('warning',['Long gap in EF ' irf_disp_iso_range(ef(iGap+[0 1],1)')])
          E3D_BASE(t_BASE(:,1)>ef(iGap,1) & t_BASE(:,1)<ef(iGap+1,1),2:4) = NaN;
          if ef(iGap+1,1)-ef(iGap,1)>DT_PC2
            irf.log('warning','splitting EF')
            T_TMP = t_BASE(t_BASE(:,1)>=ef(iGap,1) & t_BASE(:,1)<=ef(iGap+1,1));
            tint_TMP = [T_TMP(end) ttE(end,2)];
            ttE(end,2) = T_TMP(1); ttE = [ttE; tint_TMP];
          end
        end
      end
      E3D_BASE(t_BASE(:,1)<ef(1,1),2:4) = NaN; E3D_BASE(t_BASE(:,1)>ef(end,1),2:4) = NaN;

      % Construct the inertial frame
      evxb = irf_tappl(irf_cross(B_BASE,irf_resamp(V,t_BASE)),'*1e-3*(-1)');
      iE3D_BASE = E3D_BASE;
      iE3D_BASE(:,2:4) = iE3D_BASE(:,2:4) - evxb(:,2:4);
    end

    % Set Gaps to NaN. for safety we set both E and B
    if nGaps>0
      for iGap=1:nGaps
        tintGap = ttGap.TimeInterval(iGap,:);
        B_BASE(B_BASE(:,1)>=tintGap(1)-0.5 & B_BASE(:,1)<=tintGap(2)+0.5,2:end) = NaN;
        if ~isempty(ef)
          iE3D_BASE(iE3D_BASE(:,1)>=tintGap(1)-0.5 & iE3D_BASE(:,1)<=tintGap(2)+0.5,2:end) = NaN;
        end
      end
    end

    % Join intervals for BL and EF
    tt = [ttE; ttB];
    irf.log('warning',sprintf('starting with %d chunks',size(tt,1)))
    irf_disp_iso_range(tt);
    [~,ii] = sort(tt(:,1)); tt = tt(ii,:);
    idx = 1;
    while size(tt,1)>1
      ii = find(tt(:,1)>=tt(idx,1) & tt(:,1)<=tt(idx,2)); ii(ii==idx) = [];
      if isempty(ii), idx = idx+1;
        if idx>=size(tt,1), break, end
        continue,
      end
      ii = ii(1); if tt(ii,2)>tt(idx,2), tt(idx,2) = tt(ii,2); end
      tt(ii,:) = [];
    end
    irf.log('warning',sprintf('total %d chunks',size(tt,1)))
    irf_disp_iso_range(tt);

    for iChunk = 1:size(tt,1)
      tintChunk = tt(iChunk,:);
      irf.log('warning',...
        sprintf('PC12 chunk #%d : %s',iChunk,irf_disp_iso_range(tintChunk,1)))
      idxChunk = (t_BASE(:,1)>=tintChunk(1) & t_BASE(:,1)<=tintChunk(end));
      iE3D_TMP = []; if ~isempty(iE3D_BASE), iE3D_TMP = iE3D_BASE(idxChunk,:); end
      tic
      ebsp = irf_ebsp(iE3D_TMP,B_BASE(idxChunk,:),[],B0_1MIN,[],'pc12','fac',...
        'polarization','noresamp','fullB=dB','dedotb=0','nav',12,...
        'facMatrix',facMatrix);
      toc

      ebsp.r = gseR; % add position for plotting
      tlim_ebsp(tintChunk);
      if isempty(ebsp.t)
        irf.log('warning','No result for PC12'),continue
      end
      flim_ebsp(fSampB,fSampE);
      if plotFlag
        close all
        h = irf_pl_ebsp(ebsp);
        irf_zoom(h,'x',tintChunk)
        title(h(1),['THEMIS ' upper(thId) ', ' irf_disp_iso_range(tintChunk,1)])
        orient tall
        print('-dpng',['MAARBLE_TH' upper(thId) '_ULF_PC12_' irf_fname(tintChunk,5)])
      end
      if exportFlag
        maarble.export(ebsp,tintChunk,['th' thId],'pc12')
      end
    end
  end

  % Export FAC matrix
  if exportFlag
    [facMatrix.t,idxTlim]=irf_tlim(facMatrix.t,tint);
    if ~isempty(facMatrix.t)
      facMatrix.rotMatrix = facMatrix.rotMatrix(idxTlim,:,:);
      % Position is exported in GSE
      gseR_tmp = irf_resamp(gseR,facMatrix.t);
      facMatrix.r = gseR_tmp(:,2:4);
      maarble.export(facMatrix,tint,['th' thId])
    end
  end

end

  function tlim_ebsp(timeLim) % Trim ebsp to tlim
    IGNORE_FIELDS = {'f','flagFac','fullB','B0','r'};
    fieldsEBSP = fields(ebsp);
    tFields = setxor(fieldsEBSP,IGNORE_FIELDS);
    %nData = length(ebsp.t);
    [~,idx] = irf_tlim(ebsp.t,timeLim);
    for fName = tFields'
      if isempty(ebsp.(fName{:})), continue, end
      s = size(ebsp.(fName{:}));
      switch numel(s)
        case 2
          ebsp.(fName{:}) = ebsp.(fName{:})(idx,:);
        case 3
          ebsp.(fName{:}) = ebsp.(fName{:})(idx,:,:);
        otherwise
          error('wrong size!')
      end
    end
  end % tlim_ebsp()

  function flim_ebsp(fSampB,fSampE) % Thim the unmeasured frequencies
    IGNORE_FIELDS = {'t','f','flagFac','fullB','B0','r'};
    fieldsEBSP = fields(ebsp);
    tFields = setxor(fieldsEBSP,IGNORE_FIELDS);
    idxB = find(ebsp.f>fSampB/2); idxE = find(ebsp.f>fSampE/2);
    for fName = tFields'
      if isempty(ebsp.(fName{:})), continue, end
      if strcmpi(fName{:}(1:2),'ee'), idx = idxE;
      elseif strcmpi(fName{:}(1:2),'pf') %Poynting flux
        if fSampB>fSampE, idx = idxE; else, idx = idxB; end
      else, idx = idxB;
      end
      s = size(ebsp.(fName{:}));
      switch numel(s)
        case 2
          ebsp.(fName{:})(:,idx) = NaN;
        case 3
          ebsp.(fName{:})(:,idx,:) = NaN;
        otherwise
          error('wrong size!')
      end
    end
  end % flim_ebsp()
end

function data = clear_backward_jump(data,name)
if nargin <2, name = ''; else, name = [' in ' name]; end
while true
  if isempty(data), break, end
  iJump = find( diff(data(:,1))<=0 );
  if isempty(iJump), break, end
  iJump = iJump(1)+1;
  irf.log('warning',['backward time jump' name ' at ' epoch2iso(data(iJump,1))])
  ii = find(data(:,1)>=data(iJump,1)); ii(ii>=iJump) = []; data(ii,:) = [];
  irf.log('warning',sprintf('Disregarging %d points',length(ii)))
end
end
