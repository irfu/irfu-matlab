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
%  FPI IONS:
%     Vi_gse_fpi_brst, Vi_gse_fpi_brst, 'Ti_fpi_ql','Ti_fpi_brst',
%     'Ni_fpi_ql', 'Ni_fpi_brst', 'Ne_fpi_brst'
%  FGM:
%     'B_dmpa_srvy','B_gse_srvy','B_gsm_srvy',
%     'B_dmpa_brst','B_gse_brst','B_gsm_brst'
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
  'Vi_gse_fpi_ql','Ve_gse_fpi_brst','Vi_gse_fpi_brst', ...
  'Ni_fpi_ql','Ni_fpi_brst','Ne_fpi_brst',...
  'Ti_fpi_ql','Ti_fpi_brst',...
  'B_dmpa_srvy','B_gse_srvy','B_gsm_srvy','B_dmpa_brst','B_gse_brst','B_gsm_brst',...
  'dfg_ql_srvy','afg_ql_srvy'}; % XXX THESE MUST BE THE SAME VARS AS BELOW
if isempty(intersect(varStr,vars)),
  errS = ['variable not recognized: ' varStr];
  irf.log('critical',errS);
  disp('Implemented vars are:')
  for iVar = 1:length(vars)
    fprintf('  %s\n',vars{iVar})
  end
  error(errS)
end

switch varStr
  case {'R_gse','R_gsm','V_gse','V_gsm'}
    vC = varStr(1); cS = varStr(3:5);
    
    if 0 % XXX this files are WRONG!!!
    %Try to load common ephemeris file
    switch cS
      case 'gse', fSuf = '';
      case 'gsm', fSuf = 'GSM';
      otherwise, error('should not be here')
    end
    commonFile = ['/data/mms/irfu/mms' vC fSuf '.mat'];
    if exist(commonFile,'file')
      data = load(commonFile); data = data.(vC);
      data = my_tlim(data);
      if ~isempty(data)
        if mmsId==0, res = data; return; end
        res = irf.ts_vec_xyz(EpochTT(data.time),data.([cS vC mmsIdS]));
        res.time = EpochTT(res.time);
        res.name = sprintf('%s_%d',varStr,mmsId);
        res.coordinateSystem = cS;
        switch vC
          case 'R', res.units = 'km'; res.siConversion = '1e3>m';
          case 'V', res.units = 'km/s'; res.siConversion = '1e3>m/s';
          otherwise, error('should not be here')
        end
        return
      end
    end
    end
    
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
  case {'Vi_gse_fpi_ql','Vi_gse_fpi_brst','Ve_gse_fpi_brst'}
    if varStr(2)=='i', vS = 'dis';
    else vS = 'des';
    end
    if varStr(12)=='q'
      datasetName = ['mms' mmsIdS '_fpi_fast_ql_' vS];
    else
      datasetName = ['mms' mmsIdS '_fpi_brst_l1b_' vS '-moms'];
    end
    rX = mms.db_get_ts(datasetName,['mms' mmsIdS '_' vS '_bulkX'],Tint);
    if isempty(rX), return, end
    rX = comb_ts(rX);
    rY = comb_ts(mms.db_get_ts(datasetName,['mms' mmsIdS '_' vS '_bulkY'],Tint));
    rZ = comb_ts(mms.db_get_ts(datasetName,['mms' mmsIdS '_' vS '_bulkZ'],Tint));
    res = irf.ts_vec_xyz(rX.time, [rX.data rY.data rZ.data]);
    res.coordinateSystem = 'gse';
    res.name = [varStr '_' mmsIdS];
    res.units = rX.units;
    res.siConversion = rX.siConversion;
  case {'Ni_fpi_ql','Ni_fpi_brst','Ne_fpi_brst','Ti_fpi_ql','Ti_fpi_brst'}
    if varStr(2)=='i', vS = 'dis';
    else vS = 'des';
    end
    if varStr(8)=='q'
      datasetName = ['mms' mmsIdS '_fpi_fast_ql_' vS];
    else
      datasetName = ['mms' mmsIdS '_fpi_brst_l1b_' vS '-moms'];
    end
    if varStr(1)=='N' % density
      rX = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_numberDensity'],Tint);
    else % temperature
      rX = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_TempXX'],Tint);
      rY = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_TempYY'],Tint);
      rZ = mms.db_get_ts(datasetName,...
        ['mms' mmsIdS '_' vS '_TempZZ'],Tint);
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
      irf.log('notice',['Removing ',sum(ind), ...
        ' samples due to overlap AFG/DFG when transitioning between fast/slow mode.']);
      res = rTs(~ind);
    else
      res = rTs;
    end
  otherwise, error('should not be here')
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