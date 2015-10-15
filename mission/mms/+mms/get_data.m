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
%     Vi_gse_fpi_brst
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
  'Vi_gse_fpi_brst'}; % XXX THESE MUST BE THE SAME VARS AS BELOW
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
            ['mms' mmsIdS '_ql_pos_gse' cS],TintTmp);
        end
      end
      if isempty(dTmp), continue, end
      
      dTmp.data = double(dTmp.data);
      dTmpR = dTmp.resample(res.time,'spline');
      res.([cS vC mmsIdS]) = dTmpR.data; 
    end
  case 'Vi_gse_fpi_brst'
    datasetName = ['mms' mmsIdS '_fpi_brst_l1b_dis-moms'];
    rX = mms.db_get_ts(datasetName,['mms' mmsIdS '_dis_bulkX'],Tint);
    rY = mms.db_get_ts(datasetName,['mms' mmsIdS '_dis_bulkY'],Tint);
    rZ = mms.db_get_ts(datasetName,['mms' mmsIdS '_dis_bulkZ'],Tint);
    res = irf.ts_vec_xyz(rX.time, [rX.data rY.data rZ.data]);
    res.coordinateSystem = 'gse';
    res.name = [varStr '_' mmsIdS];
    res.units = rX.units;
    res.siConversion = rX.siConversion;
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