function res = get_data(varStr, Tint, mmsId)
%MMS.GET_DATA  Load a variable
%
%  res = MMS.GET_DATA(varStr, Tint, [mmsId])
%
%  mmsId=0 (DEFAULT) means 1:4
%
%  varStr is one of:
%  EPHEMERIS:
%     R_gse, R_gsm, V_gse, V_gsm
%     tetra_quality
%  FPI IONS:
%     'Ni_fpi_brst_l2' (alias:'Ni_fpi_brst'), 'Ni_fpi_fast_l2',...
%     'Ni_fpi_sitl','Ni_fpi_ql',...
%     'Ni_fpi_brst_l1b','Ni_fpi_fast_l1b',...
%     'Vi_dbcs_fpi_brst_l2' (alias:'Vi_dbcs_fpi_brst'), 'Vi_dbcs_fpi_fast_l2',...
%     'Vi_gse_fpi_sitl', 'Vi_gse_fpi_ql',...
%     'Vi_gse_fpi_brst_l1b','Vi_gse_fpi_fast_l1b',...
%     'Tsi_fpi_brst_l2' (alias:'Tsi_fpi_brst'), 'Tsi_fpi_fast_l2',... %scalar temperature
%     'Tsi_fpi_sitl','Tsi_fpi_ql',...
%     'Tsi_fpi_brst_l1b','Tsi_fpi_fast_l1b',...
%     'Ti_fpi_ql','Ti_fpi_brst', 'Ti_fpi_fast_l2'
%  FPI ELECTRONS:
%     'Ne_fpi_brst_l2' (alias:'Ne_fpi_brst'),'Ne_fpi_fast_l2',...
%     'Ne_fpi_sitl','Ne_fpi_ql',...
%     'Ne_fpi_brst_l1b','Ne_fpi_fast_l1b',...
%     'Ve_dbcs_fpi_brst_l2' (alias:'Ve_dbcs_fpi_brst'), 'Ve_dbcs_fpi_fast_l2',...
%     'Ve_gse_fpi_sitl', 'Ve_gse_fpi_ql',...
%     'Ve_gse_fpi_brst_l1b','Ve_gse_fpi_fast_l1b',...
%     'Tse_fpi_brst_l2' (alias:'Tse_fpi_brst'), 'Tse_fpi_fast_l2',... %scalar temperature
%     'Tse_fpi_sitl','Tse_fpi_ql',...
%     'Tse_fpi_brst_l1b','Tse_fpi_fast_l1b',...
%     Loads into tensor of order 2:
%     'Te_fpi_ql','Te_fpi_brst','Te_fpi_brst_l2', 'Te_fpi_fast_l2'
%     'Pe_fpi_ql','Pe_fpi_brst','Pe_fpi_brst_l2'
%  FGM:
%     'B_dmpa_srvy','B_gse_srvy','B_gsm_srvy',
%     'B_dmpa_brst','B_gse_brst','B_gsm_brst'
% Added 2016-02-18 by Per-Arne:
%  HPCA:
%     'Nhplus_hpca_sitl'
%
% Example:
%   Tint = irf.tint('2015-09-21T00:00:00Z/2015-09-21T17:00:00Z');
%   V1 = mms.get_data('V_gse',Tint,1); % SC GSE velocity for MMS1
%   R  = mms.get_data('R_gse',Tint);   % SC GSE position for all MMS SC

res = [];

if nargin<3, mmsId = 0; end
if isempty(intersect(mmsId,0:4)),
  errS = ['invalid MMS ID: ' mmsId];
  irf.log('critical',errS); error(errS)
end
mmsIdS = num2str(mmsId);

if ~isa(Tint,'GenericTimeArray')
  errS = 'TINT must be of GenericTimeArray type';
  irf.log('critical',errS); error(errS)
elseif Tint.stop-Tint.start<=0,
  errS = 'TINT duration is zero or negative';
  irf.log('critical',errS); error(errS)
end

vars = {'R_gse','R_gsm','V_gse','V_gsm',...
  'Vi_dbcs_fpi_brst_l2', 'Vi_dbcs_fpi_brst', 'Vi_dbcs_fpi_fast_l2',...
  'Vi_gse_fpi_sitl', 'Vi_gse_fpi_ql',...
  'Vi_gse_fpi_brst_l1b','Vi_gse_fpi_fast_l1b',...
  'Ve_dbcs_fpi_brst_l2','Ve_dbcs_fpi_brst', 'Ve_dbcs_fpi_fast_l2',...
  'Ve_gse_fpi_sitl', 'Ve_gse_fpi_ql',...
  'Ve_gse_fpi_brst_l1b','Ve_gse_fpi_fast_l1b',...
  'Ni_fpi_brst_l2','Ni_fpi_brst','Ni_fpi_fast_l2',...
  'Ni_fpi_sitl','Ni_fpi_ql',...
  'Ni_fpi_brst_l1b','Ni_fpi_fast_l1b',...
  'Ne_fpi_brst_l2','Ne_fpi_brst','Ne_fpi_fast_l2',...
  'Ne_fpi_sitl','Ne_fpi_ql',...
  'Ne_fpi_brst_l1b','Ne_fpi_fast_l1b',...
  'Pe_fpi_ql','Pe_fpi_brst','Pe_fpi_brst_l2','Pi_fpi_brst_l2',...
  'Tsi_fpi_brst_l2','Tsi_fpi_brst','Tsi_fpi_fast_l2',...
  'Tsi_fpi_sitl','Tsi_fpi_ql',...
  'Tsi_fpi_brst_l1b','Tsi_fpi_fast_l1b',...
  'Tse_fpi_brst_l2','Tse_fpi_brst','Tse_fpi_fast_l2',...
  'Tse_fpi_sitl','Tse_fpi_ql',...
  'Tse_fpi_brst_l1b','Tse_fpi_fast_l1b',...
  'Ti_dbcs_fpi_brst_l2','Ti_dbcs_fpi_brst','Ti_dbcs_fpi_fast_l2',...
  'Ti_gse_fpi_sitl','Ti_gse_fpi_ql',...
  'Ti_gse_fpi_brst_l1b','Ti_gse_fpi_fast_l1b',...
  'Te_dbcs_fpi_brst_l2','Te_dbcs_fpi_brst','Te_dbcs_fpi_fast_l2',...
  'Te_gse_fpi_sitl','Te_gse_fpi_ql',...
  'Te_gse_fpi_brst_l1b','Te_gse_fpi_fast_l1b',...
  'Pi_dbcs_fpi_brst_l2','Pi_dbcs_fpi_brst','Pi_dbcs_fpi_fast_l2',...
  'Pi_gse_fpi_sitl','Pi_gse_fpi_ql',...
  'Pi_gse_fpi_brst_l1b','Pi_gse_fpi_fast_l1b',...
  'Pe_dbcs_fpi_brst_l2','Pe_dbcs_fpi_brst','Pe_dbcs_fpi_fast_l2',...
  'Pe_gse_fpi_sitl','Pe_gse_fpi_ql',...
  'Pe_gse_fpi_brst_l1b','Pe_gse_fpi_fast_l1b',...
  'B_dmpa_srvy','B_gse_srvy','B_gsm_srvy','B_dmpa_brst','B_gse_brst','B_gsm_brst',...
  'dfg_ql_srvy','afg_ql_srvy','tetra_quality',...
  'Nhplus_hpca_sitl'}; % XXX THESE MUST BE THE SAME VARS AS BELOW
if isempty(intersect(varStr,vars)),
  errS = ['variable not recognized: ' varStr];
  irf.log('critical',errS);
  vars = sort(vars);
  disp('Implemented vars are:')
  for iVar = 1:length(vars)
    fprintf('  %s\n',vars{iVar})
  end
  error(errS)
end

switch varStr
  case {'R_gse','R_gsm','V_gse','V_gsm'}
    vC = varStr(1); cS = varStr(3:5);
    
    if mmsId>0
      res = mms.db_get_ts(['mms' mmsIdS '_mec_srvy_l2_epht89d'],...
        ['mms' mmsIdS '_mec_' lower(vC) '_' cS],Tint);
      if ~isempty(res), return, end
      
      % LAST RESORT: Load position of MAG files
      if vC=='R'
      % Load from L2pre B
      res = mms.db_get_ts(...
        ['mms' mmsIdS '_dfg_srvy_l2pre'],['mms' mmsIdS '_pos_' cS],Tint);
      if ~isempty(res), return, end
      % Load from QL B
      res = mms.db_get_ts(...
        ['mms' mmsIdS '_dfg_srvy_ql'],['mms' mmsIdS '_ql_pos_' cS],Tint);
      return
      end
    end
    
    % Do resampling similar to mms_update_ephemeris
    TintTmp = irf.tint(Tint.start+(-60),Tint.stop+60);
    TintTmp = EpochUnix([fix(TintTmp.start.epochUnix/60)*60 ...
      ceil(TintTmp.stop.epochUnix/60)*60]);
    TintTmp = EpochTT(TintTmp);
    res.time = EpochTT((TintTmp.start.epoch:int64(30*1e9):TintTmp.stop.epoch)');
    for mmsId=1:4
      mmsIdS = num2str(mmsId);
      dTmp = mms.db_get_ts(['mms' mmsIdS '_mec_srvy_l2_epht89d'],...
        ['mms' mmsIdS '_mec_' lower(vC) '_' cS],Tint);
      if isempty(dTmp) &&  vC=='V', continue, end
      
      if isempty(dTmp)
        % LAST RESORT: Load position of MAG files
        % Load from L2pre B
        dTmp = mms.db_get_ts(['mms' mmsIdS '_dfg_srvy_l2pre'],...
          ['mms' mmsIdS '_pos_' cS],TintTmp);
        if isempty(dTmp)
          % Load from QL B
          dTmp = mms.db_get_ts(['mms' mmsIdS '_dfg_srvy_ql'],...
            ['mms' mmsIdS '_ql_pos_' cS],TintTmp);
        end
      end
      if isempty(dTmp), continue, end
      dTmp = comb_ts(dTmp);
      dTmp.data = double(dTmp.data);
      dTmpR = dTmp.resample(res.time,'spline');
      res.([cS vC mmsIdS]) = dTmpR.data; 
    end 
end

Vr = splitVs(varStr);
datasetName = ['mms' mmsIdS '_' Vr.inst '_' Vr.tmmode '_' Vr.lev];
compS = ''; pref = ''; suf = '';

switch Vr.inst
  case 'fpi'
    switch Vr.param(end)
      case 'i', sensor = 'dis';
      case 'e', sensor = 'des';
      otherwise 
        error('invalid specie')
    end
    
    switch Vr.lev
      case {'l2','l2pre','l1b'}
        datasetName = [datasetName '_' sensor '-moms'];
      case 'ql'
        datasetName = [datasetName '_' sensor];
      case 'sitl'
      otherwise, error('should not be here')
    end
    
    switch Vr.param
      case {'Ni','Ne'}
        switch Vr.lev
          case {'l2','l2pre'}
            pref = ['mms' mmsIdS '_' sensor '_numberdensity_dbcs_' Vr.tmmode];
          case 'l1b'
            pref = ['mms' mmsIdS '_' sensor '_numberdensity'];
          case 'ql'
            pref = ['mms' mmsIdS '_' sensor '_numberDensity'];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'numberDensity'];
          otherwise, error('should not be here')
        end
        res = get_ts('scalar');
        return
      case {'Tsi','Tse'}
        getQ = 'trace';
        switch Vr.lev
          case {'l2','l2pre'}
            pref = ['mms' mmsIdS '_' sensor '_temp'];
            suf = ['_dbcs_' Vr.tmmode];
            compS = struct('xx','xx','yy','yy','zz','zz');
          case {'l1b','ql'}
            pref = ['mms' mmsIdS '_' sensor '_Temp'];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'temp'];
            getQ = 'ts';
          otherwise, error('should not be here')
        end
        res = get_ts(getQ);
        return
      case {'Ti', 'Te', 'Pi', 'Pe'}
        switch Vr.param(1)
          case 'T' % temperature
            momType = 'Temp';
          case 'P' % pressure
            momType = 'Pres';
          otherwise, error('should not be here 2')
        end
        switch Vr.lev
          case {'l2','l2pre'}
            pref = ['mms' mmsIdS '_' sensor '_' lower(momType)];
            suf = ['_' Vr.cs '_' Vr.tmmode];
            compS = struct('xx','xx','xy','xy','xz','xz','yy','yy','yz','yz','zz','zz'); 
          case {'l1b','ql'}
            pref = ['mms' mmsIdS '_' sensor '_' momType];
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' upper(sensor) 'temp'];
            getQ = 'ts';
          otherwise, error('should not be here')
        end
        res = get_ts('tensor2');
        return
      case {'Vi','Ve'}
        pref = ['mms' mmsIdS '_' sensor '_bulk'];
        switch Vr.lev
          case {'l2','l2pre'}
            suf = ['_' Vr.cs '_' Vr.tmmode];
            compS = struct('x','x','y','y','z','z');
          case 'l1b'
          case 'ql'
          case 'sitl'
            pref = ['mms' mmsIdS '_fpi_' Vr.param(end) 'BulkV_'];
            suf = '_DSC';
          otherwise, error('should not be here')
        end
        res = get_ts('vector');
        return
    end    
  otherwise
    %error('not implemented yet')
end

switch varStr
  case {'Pe_fpi_ql','Pe_fpi_brst','Te_fpi_ql',... % QL
      'Te_fpi_brst','Te_fpi_fast',... % L1b
      'Pe_fpi_brst_l2','Te_fpi_brst_l2','Pi_fpi_brst_l2','Ti_fpi_brst_l2',...
      'Pe_fpi_fast_l2','Te_fpi_fast_l2','Pi_fpi_fast_l2','Ti_fpi_fast_l2'
      }
    isL2 = 0;
    if varStr(2)=='i', vS = 'dis';
    else vS = 'des';
    end
    if varStr(8)=='q'
      datasetName = ['mms' mmsIdS '_fpi_fast_ql_' vS];
    else
      if varStr(8)=='b', modeS = 'brst'; else modeS = 'fast'; end
      if varStr(end) == '2'
        isL2 = 1;
        datasetName = ['mms' mmsIdS '_fpi_' modeS '_l2_' vS '-moms'];
      else
        datasetName = ['mms' mmsIdS '_fpi_' modeS '_l1b_' vS '-moms'];
      end
    end
    if varStr(1)=='T' % temperature
      momType = 'Temp';
    elseif varStr(1)=='P' % pressure
      momType = 'Pres';
    else return;
    end
    
    if isL2 % L2 data has different variable names
      rXX = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('XX') '_dbcs_' modeS],Tint);
      if isempty(rXX), return, end
      rXY = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('XY') '_dbcs_' modeS],Tint);
      rXZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('XZ') '_dbcs_' modeS],Tint);
      rYY = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('YY') '_dbcs_' modeS],Tint);
      rYZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('YZ') '_dbcs_' modeS],Tint);
      rZZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' lower(momType) lower('ZZ') '_dbcs_' modeS],Tint); 
      coordinateSystem = 'DBCS';
    else
      rXX = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'XX'],Tint);
      if isempty(rXX), return, end
      rXY = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'XY'],Tint);
      rXZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'XZ'],Tint);
      rYY = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'YY'],Tint);
      rYZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'YZ'],Tint);
      rZZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_' momType 'ZZ'],Tint);  
      coordinateSystem = '';
    end 
    
    rData = nan(rXX.length,3,3);
    rData(:,1,1) = rXX.data;
    rData(:,1,2) = rXY.data;
    rData(:,1,3) = rXZ.data;
    rData(:,2,1) = rXY.data;
    rData(:,2,2) = rYY.data;
    rData(:,2,3) = rYZ.data;
    rData(:,3,1) = rXZ.data;
    rData(:,3,2) = rYZ.data;
    rData(:,3,3) = rZZ.data;
        
    res = irf.ts_tensor_xyz(rXX.time, rData);
    res.name = [varStr '_' mmsIdS];
    res.units = rXX.units;
    res.siConversion = rXX.siConversion;
    res.coordinateSystem = coordinateSystem;
  case {'Ti_fpi_ql','Ti_fpi_brst'}
    if varStr(2)=='i', vS = 'dis';
    else vS = 'des';
    end
    switch varStr(8)
      case 's', datasetName = ['mms' mmsIdS '_fpi_fast_sitl']; vS = 'fpi_DIS';
      case 'q', datasetName = ['mms' mmsIdS '_fpi_fast_ql_' vS]; vS = [vS '_'];
      case 'f'
      otherwise
        datasetName = ['mms' mmsIdS '_fpi_brst_l1b_' vS '-moms']; vS = [vS '_'];
    end
    pref = ['mms' mmsIdS '_' vS];
    if varStr(1)=='N' % density
      rX = mms.db_get_ts(datasetName,[pref 'numberDensity'],Tint);
      rX = comb_ts(rX);
    else % temperature
      rX = mms.db_get_ts(datasetName, [pref 'TempXX'],Tint); rX = comb_ts(rX);
      rY = mms.db_get_ts(datasetName, [pref 'TempYY'],Tint); rY = comb_ts(rY);
      rZ = mms.db_get_ts(datasetName, [pref 'TempZZ'],Tint); rZ = comb_ts(rZ);
      rX.data = rX.data + rY.data + rZ.data;
    end
    if isempty(rX), return, end
    res = irf.ts_scalar(rX.time, rX.data);
    res.name = [varStr '_' mmsIdS];
    res.units = rX.units;
    res.siConversion = rX.siConversion;
  case {'Ni_fpi_fast_l2','Ne_fpi_fast_l2'}
    if varStr(2)=='i', vS = 'dis';
    else vS = 'des';
    end
    datasetName = ['mms' mmsIdS '_fpi_fast_l2_' vS '-moms']; vS = [vS '_'];
    pref = ['mms' mmsIdS '_' vS];
    if varStr(1)=='N' % density
      rX = mms.db_get_ts(datasetName,[pref 'numberdensity_dbcs_fast'],Tint);
      rX = comb_ts(rX);
    else % temperature
      rX = mms.db_get_ts(datasetName, [pref 'tempxx_dbcs_fast'],Tint); rX = comb_ts(rX);
      rY = mms.db_get_ts(datasetName, [pref 'tempyy_dbcs_fast'],Tint); rY = comb_ts(rY);
      rZ = mms.db_get_ts(datasetName, [pref 'tempzz_dbcs_fast'],Tint); rZ = comb_ts(rZ);
      rX.data = rX.data + rY.data + rZ.data;
    end
    if isempty(rX), return, end
    res = irf.ts_scalar(rX.time, rX.data);
    res.name = [varStr '_' mmsIdS];
    res.units = rX.units;
    res.siConversion = rX.siConversion;
  case {'B_dmpa_srvy','B_gse_srvy','B_gsm_srvy','B_dmpa_brst','B_gse_brst','B_gsm_brst'}
    instr = 'dfg';
    tk = tokenize(varStr,'_'); cS = tk{2}; dLev = tk{3};
    datasetName = ['mms' mmsIdS '_' instr '_' dLev '_l2pre'];
    varName = ['mms' mmsIdS '_' instr '_' dLev '_l2pre_' cS];
    rTs = mms.db_get_ts(datasetName, varName, Tint);
    if isempty(rTs), return, end
    if strcmp(dLev,'srvy')
      ind = diff(rTs.time.ttns) <= 122000; % FIXME: what is brst min dt for A/DFG?
      if( sum(ind) < (length(rTs)-2) )
        % Remove samples that are too close, but ensure some output if only
        % two samples with very high sample rate.
        irf.log('notice',['Removing ',sum(ind), ...
          ' samples due to overlap AFG/DFG when transitioning between fast/slow mode.']);
        res = rTs(~ind);
      else res = rTs;
      end
    end
  case {'dfg_ql_srvy', 'afg_ql_srvy'} % FIXME: Correct name, and above as well!!
    instr = varStr(1:3);
    datasetName = ['mms', mmsIdS, '_', instr, '_srvy_ql'];
    varName = ['mms', mmsIdS, '_', instr, '_srvy_dmpa'];
    rTs = mms.db_get_ts(datasetName, varName, Tint);
    if isempty(rTs), return, end
    rTs = comb_ts(rTs); % Try to combine multiple results.
    ind = diff(rTs.time.ttns) <= 122000; % FIXME: what is brst min dt for A/DFG?
    if( sum(ind) < (length(rTs)-2) )
      % Remove samples that are too close, but ensure some output if only
      % two samples with very high sample rate.
      irf.log('notice',['Removing ',num2str(sum(ind)), ...
        ' samples due to overlap AFG/DFG when transitioning between fast/slow mode.']);
      res = rTs(~ind);
    else
      res = rTs;
    end
  case 'tetra_quality'
    % Begin looking for Def. quality
    quality = mms.db_get_variable('mms_ancillary_defq','quality',Tint);
    if isempty(quality)
      irf.log('warning', 'Did not find any definite tetrahedra quality. Looking for predicted.');
      list = mms.db_list_files('mms_ancillary_predq',Tint);
      if(~isempty(list))
        % Load the last predicted file to match Tint
        quality = mms_load_ancillary([list(end).path, filesep, list(end).name], 'predq');
      end
    end
    if(~isempty(quality))
      rTs = irf.ts_scalar(EpochTT(quality.time), quality.quality);
      res = rTs.tlim(Tint);
    end
  case {'Nhplus_hpca_sitl'}
    datasetName = ['mms' mmsIdS '_hpca_srvy_sitl_moments'];
    rX = mms.db_get_ts(datasetName,...
      ['mms' mmsIdS '_hpca_hplus_number_density'],Tint);
    if isempty(rX), return, end
    res = irf.ts_scalar(rX.time, rX.data);
    res.name = [varStr '_' mmsIdS];
    res.units = rX.units;
    res.siConversion = rX.siConversion;
  otherwise, error('should not be here')
end

  function res = get_ts(dataType)
    res = [];
    switch dataType
      case 'scalar'
        rX = mms.db_get_ts(datasetName,pref,Tint);
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' datasetName '(' pref ')'])
          return
        end
        rX = comb_ts(rX);
        res = irf.ts_scalar(rX.time, rX.data);
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'vector'
        if isempty(compS), compS.x = 'X'; compS.y = 'Y'; compS.z = 'Z'; end
        rX = mms.db_get_ts(datasetName,[pref compS.x suf],Tint);
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' datasetName '(' [pref compS.x suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = comb_ts(mms.db_get_ts(datasetName,[pref compS.y suf],Tint));
        rZ = comb_ts(mms.db_get_ts(datasetName,[pref compS.z suf],Tint));
        res = irf.ts_vec_xyz(rX.time, [rX.data rY.data rZ.data]);
        res.coordinateSystem = Vr.cs;
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'trace'
        if isempty(compS), compS.xx = 'XX'; compS.yy = 'YY'; compS.zz = 'ZZ'; end
        rX = mms.db_get_ts(datasetName, [pref compS.xx suf],Tint); 
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' datasetName '(' [pref compS.xx suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = mms.db_get_ts(datasetName, [pref compS.yy suf],Tint); rY = comb_ts(rY);
        rZ = mms.db_get_ts(datasetName, [pref compS.yy suf],Tint); rZ = comb_ts(rZ);
        rX.data = rX.data + rY.data + rZ.data;
        res = irf.ts_scalar(rX.time, rX.data);
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'ts'
        if isempty(compS), compS.par = 'Para'; compS.perp = 'Perp'; end
        rX = mms.db_get_ts(datasetName, [pref compS.par suf],Tint); 
        if isempty(rX)
          irf.log('warning',...
            ['No data for ' datasetName '(' [pref compS.par suf] ')'])
          return
        end
        rX = comb_ts(rX);
        rY = mms.db_get_ts(datasetName, [pref compS.perp suf],Tint); rY = comb_ts(rY);
        rX.data = rX.data/3 + rY.data*2/3;
        res = irf.ts_scalar(rX.time, rX.data);
        res.name = [varStr '_' mmsIdS];
        res.units = rX.units;
        res.siConversion = rX.siConversion;
      case 'tensor2'
        if isempty(compS), 
          compS = struct('xx','XX','xy','XY','xz','XZ','yy','YY','yz','YZ','zz','ZZ'); 
        end
        rXX = mms.db_get_ts(datasetName,[pref compS.xx suf],Tint);
        if isempty(rXX),irf.log('warning',...
            ['No data for ' datasetName '(' [pref compS.par suf] ')'])
          return
        end
        rXX = comb_ts(rXX);
        rXY = mms.db_get_ts(datasetName,[pref compS.xy suf],Tint);rXY = comb_ts(rXY);
        rXZ = mms.db_get_ts(datasetName,[pref compS.xz suf],Tint);rXZ = comb_ts(rXZ);
        rYY = mms.db_get_ts(datasetName,[pref compS.yy suf],Tint);rYY = comb_ts(rYY);
        rYZ = mms.db_get_ts(datasetName,[pref compS.yz suf],Tint);rYZ = comb_ts(rYZ);
        rZZ = mms.db_get_ts(datasetName,[pref compS.zz suf],Tint);rZZ = comb_ts(rZZ);
    
        rData = nan(rXX.length,3,3);
        rData(:,1,1) = rXX.data;
        rData(:,1,2) = rXY.data;
        rData(:,1,3) = rXZ.data;
        rData(:,2,1) = rXY.data;
        rData(:,2,2) = rYY.data;
        rData(:,2,3) = rYZ.data;
        rData(:,3,1) = rXZ.data;
        rData(:,3,2) = rYZ.data;
        rData(:,3,3) = rZZ.data;
        
        res = irf.ts_tensor_xyz(rXX.time, rData);
        res.name = [varStr '_' mmsIdS];
        res.units = rXX.units;
        res.siConversion = rXX.siConversion;
        res.coordinateSystem = Vr.cs;
      otherwise
        error('data type not implemented')
    end
  end
  function d = my_tlim(d)
    idx = tlim(EpochTT(d.time),Tint);
    f = fields(d);
    for iF = 1:length(f)
      d.(f{iF}) = d.(f{iF})(idx,:);
    end
  end
end %% MAIN

function TsOut = comb_ts(TsIn)
 if ~iscell(TsIn), TsOut = TsIn; return; end
 TsOut = TsIn{1};
 for i=2:numel(TsIn)
  TsOut = combine(TsOut, TsIn{i});
 end
end

function Res = splitVs(varStr)

tk = tokenize(varStr,'_');
nTk = length(tk);
if nTk <3 || nTk > 5, error('invalig STRING format'), end

param = tk{1};
switch param
  case {'Ni', 'Ne', 'Nhplus', 'Tsi', 'Tse'}
    tensorOrder = 0;
  case {'Vi', 'Ve', 'B', 'E'}
    tensorOrder = 1;
  case {'Pi', 'Pe', 'Ti', 'Te'}
    tensorOrder = 2;
  otherwise 
    error('invalid PARAM')
end

coordinateSystem = []; idx = 1;
if tensorOrder > 0
  coordinateSystem = tk{idx+1}; idx = idx + 1;
  switch coordinateSystem
    case {'gse','gsm','dsl','dbcs','dmpa'}
    otherwise
      error('invalid COORDINATE_SYS')
  end
end

instrument = tk{idx+1}; idx = idx + 1;
switch instrument
  case {'fpi','edp','edi','hpca','fgm','scm'}
  otherwise
    switch param
      case B, instrument = 'fgm'; idx = idx - 1;
      case {'Ve','Te','Ne','Pe'}, instrument = 'fpi'; idx = idx - 1;
      otherwise
        error('invalid INSTRUMENT')
    end
end

tmMode = tk{idx+1}; idx = idx + 1;
switch tmMode
  case {'brst','fast','slow','srvy'}
  otherwise
    tmMode = 'fast'; idx = idx - 1;
    irf.log('warning','assuming TM_MODE = FAST')
end

if length(tk)==idx, dataLevel = 'l2'; %default
else
  dataLevel = tk{idx+1};
  switch dataLevel
    case {'ql','sitl','l1b','l2a','l2pre','l2'}
    otherwise
      error('invalid DATA_LEVEL level')
  end
end

Res = struct('param',param,'to',tensorOrder,'cs',coordinateSystem,...
  'inst',instrument,'tmmode',tmMode,'lev',dataLevel);
end
