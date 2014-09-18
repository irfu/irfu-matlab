function [ VarOut ] = mms_sdc_sdp_datamanager( param, dataObj )
% mms_sdc_sdp_datamanager will store and retrive data for
%  	mms_sdc_sdp_datamanager( dataType, dataObj ) will store
%  	appropriate data variables related to dataType in the global variable
%  	DATAC for later retrival.
%  
%   [varOut] = mms_sdc_sdp_datamanager( variable ) will return the variable
%   requested to varOut, if no such variable has been created already it
%   will try and calculate it using the stored data.
%   
%  	Example:
%
%  		mms_sdc_sdp_datamanager('dce',dceDataObj)
%  		mms_sdc_sdp_datamanager('phase')
%  
%   	See also DATAOBJ, MMS_SDC_SDP_CDF_IN_PROCESS.

narginchk(1,2); % One argument is simply retreive variable, Two arguments
% store "dataObj" as "dataType".

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
global DATAC; % Here we store read data.

if ~ischar(param),
  errStr = 'PARAM must be a string';
  irf.log('critical', errStr); error(err_str);
end

if strcmpi(param, 'init')
  % Initialize
  if nargin==1
    errStr = 'INIT requires second input argument';
    irf.log('critical', errStr); error(errStr);
  elseif ~isstruct(dataObj)
    errStr = 'Second input argument for INIT must be a structure';
    irf.log('critical', errStr); error(errStr);
  elseif ~isfield(dataObj,'scId') || ~isnumeric(dataObj.scId) || ...
      isempty(intersect(dataObj.scId, MMS_CONST.MMSids))
    errStr = 'Invalid input for init_struct.scId';
    irf.log('critical', errStr); error(errStr);
  end
  DATAC.scId = dataObj.scId;
  if ~isfield(dataObj,'tmMode')
    DATAC.tmMode = 1;
    irf.log('warining',['init_struct.tmMode not specified, defaulting to '''...
      MMS_CONST.TmModes{DATAC.tmMode} ''''])
  elseif ~isnumeric(dataObj.tmMode) || ...
      isempty(intersect(dataObj.tmMode, 1:numel(MMS_CONST.TmModes)))
    errStr = 'Invalid input for init_struct.tmMode';
    irf.log('critical', errStr); error(errStr);
  else DATAC.tmMode = dataObj.tmMode;
  end
  if ~isfield(dataObj,'procId')
    DATAC.procId = 1;
    irf.log('warining',['init_struct.procId not specified, defaulting to '''...
      MMS_CONST.SDCProcs{DATAC.procId} ''''])
  elseif ~isnumeric(dataObj.procId) || ...
      isempty(intersect(dataObj.procId, 1:numel(MMS_CONST.SDCProcs)))
    errStr = 'Invalid input for init_struct.tmMode';
    irf.log('critical', errStr); error(errStr);
  else DATAC.procId = dataObj.procId;
  end
  DATAC.dce = [];
  DATAC.dcv = [];
  DATAC.hk_101 = [];
  return
end

if ~isfield(DATAC, 'scId')
  err_str = 'Data mamager not initialized! Run: mms_sdc_sdp_datamanager(''init'',init_struct)';
  irf.log('critical', err_str);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
end

param = lower(param);

if(nargin==2)   
    % Make sure first argument is a dataobj class object, 
    % otherwise a read cdf file.
    if isa(dataObj,'dataobj') % do nothing
    elseif ischar(dataObj) && exist(dataObj, 'file')
        % If it is not a read cdf file, is it an unread cdf file? Read it.
        irf.log('warning',['Loading ' param ' from file: ', dataObj]);
        dataObj = dataobj(dataObj, 'KeepTT2000');
    else
        err_str = 'MMS_SDC_SDP_DATAMANAGER unknown input arguments.';
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);  
    end
    
    if( isfield(DATAC, param) ) && ~isempty(DATAC.(param))
      % Error, Warning or Notice for replacing the data variable?
      err_str = ['MMS_SDC_SDP_DATAMANAGER replacing existing ', ...
        'variable with new data in ' param];
      irf.log('critical', err_str);
      error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    
    varPrefix = sprintf('mms%d_sdp_',DATAC.scId);
    
    switch(param)
      case('dce')
        sensors = {'e12','e34','e56'};
        init_param()
        
      case('dcv')
        sensors = {'v1','v2','v3','v4','v5','v6'};
        init_param()
        v_from_e_and_v()
        
      case('hk_101')
        varPrefix = sprintf('mms%d_101_',DATAC.scId);
        DATAC.(param) = [];
        DATAC.(param).dataObj = dataObj;
        x = getdep(dataObj,[varPrefix 'cmdexec']);
        DATAC.(param).time = x.DEPEND_O.data;
        check_monoton_timeincrease(DATAC.(param).time, param);
        % Add sunpulse times (TT2000) of last recieved sunpulse.
        DATAC.(param).sunpulse = dataObj.data.([varPrefix 'sunpulse']).data;
        % Add sunpulse indicator, real: 0, SC pseudo: 1, CIDP pseudo: 2.
        DATAC.(param).sunssps = dataObj.data.([varPrefix 'sunssps']).data;
        % Add CIDP sun period (in microseconds, 0 if sun pulse not real.
        DATAC.(param).iifsunper = dataObj.data.([varPrefix 'iifsunper']).data;
      otherwise
        % Not yet implemented.
        err_str = ['MMS_SDC_SDP_DATAMANAGER unknown second ', ...
          'parameter. The datatype ', param, ...
          ' is not yet implemented'];
        irf.log('critical',err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
        
    end
    
elseif nargin==1
        switch( param )
            case('dcephase')
                % Phase, from sunpulse for now.
                if isfield(DATAC,'dce') && ...
                    isfield(DATAC.dce,'phase') && ...
                    ~isempty(DATAC.dce.phase)
                  VarOut = DATAC.dce.phase;
                else
                    % Calculate it, store it and return variable.
                    DATAC.dce.phase = mms_sdc_sdp_phase( ...
                        DATAC.dce.time, ...
                        DATAC.hk_101.sunpulse);
                    VarOut = DATAC.dce.phase;
                end
                
                
                %%
                %THIS SHOULD BE CHANGED BUT FOR TESTING AND VERIFICATION
                %PURPOSES IT IS DONE THIS WAY. THE FOLLOWING MUST BE
                %CHANGED.
            case('dcetime')
                % Timestamp of dce
                if( isfield(DATAC.dce, 'time') && ...
                        ~isempty(DATAC.dce.time) )
                    VarOut = DATAC.dce.time;
                else
                    % Error
                    err_str = ['The requested variable ', param, ...
                        'does not exist.'];
                    irf.log('critical',err_str);
                    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
                end
                
            case('dcvtime')
                % Timestamp of dce
                if( isfield(DATAC.dcv, 'time') && ...
                        ~isempty(DATAC.dcv.time) )
                    VarOut = DATAC.dcv.time;
                else
                    % Error
                    err_str = ['The requested variable ', param, ...
                        'does not exist.'];
                    irf.log('critical',err_str);
                    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
                end
                
                %% END OF TESTING AND VERIFICATION CODE TO BE CHANGED
                
                
            otherwise
                % FIXME: Not yet implemented.
                err_str = 'MMS_SDC_SDP_DATAMANAGER variable not yet implemented.';
                irf.log('critical', err_str);
                error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
        end
end

  function v_from_e_and_v
    % Compute V from E and the other V
    % typical situation is V2 off, V1 on
    % E12[mV/m] = ( V1[V] - V2[V] ) / L[km]
    NOM_BOOM_L = .12; % 120 m
    MSK_OFF = MMS_CONST.Bitmask.SIGNAL_OFF;
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen}; senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      senA_off = bitand(DATAC.dcv.(senA).bitmask, MSK_OFF);
      senB_off = bitand(DATAC.dcv.(senB).bitmask, MSK_OFF);
      idxOneSig = xor(senA_off,senB_off);
      if ~any(idxOneSig), return, end
      iVA = idxOneSig & ~senA_off;
      if any(iVA),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senB,senA,senE,sum(iVA)))
        DATAC.dcv.(senB).data(iVA) = DATAC.dcv.(senA).data(iVA) - ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVA);
      end
      iVB = idxOneSig & ~senB_off;
      if any(iVB),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senA,senB,senE,sum(iVA)))
        DATAC.dcv.(senA).data(iVB) = DATAC.dcv.(senB).data(iVB) + ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVB);
      end
    end
  end

  function init_param
    DATAC.(param) = [];
    if ~all(diff(dataObj.data.([varPrefix 'samplerate_' param]).data)==0)
      err_str = 'MMS_SDC_SDP_DATAMANAGER changing sampling rate not yet implemented.';
      irf.log('warning', err_str);
      %error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    DATAC.(param).dataObj = dataObj;
    fileVersion = DATAC.(param).dataObj.GlobalAttributes.Data_version{:};
    DATAC.(param).fileVersion = struct(...
      'major', str2double(fileVersion(2)),...
      'minor', str2double(fileVersion(4)),...
      'revision', str2double(fileVersion(4)));
    x = getdep(dataObj,[varPrefix param '_sensor']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    sensorData = dataObj.data.([varPrefix param '_sensor']).data;
    if isempty(sensors), return, end
    probeEnabled = resample_probe_enable(sensors);
    for iSen=1:numel(sensors)
      DATAC.(param).(sensors{iSen}) = struct(...
        'data',sensorData(:,iSen), ...
        'bitmask',zeros(size(sensorData(:,iSen))));
      %Set disabled bit
      idxDisabled = probeEnabled(:,iSen)==0;
      DATAC.(param).(sensors{iSen}).bitmask(idxDisabled) = ...
        bitor(DATAC.(param).(sensors{iSen}).bitmask(idxDisabled), ...
        MMS_CONST.Bitmask.SIGNAL_OFF);
      DATAC.(param).(sensors{iSen}).data(idxDisabled,:) = NaN;
    end
  end

  function res = resample_probe_enable(fields)
  % resample probe_enabled data to E-field cadense
    probe = fields{1};
    flag = get_variable(dataObj,[varPrefix probe '_enable']);
    dtSampling = median(diff(flag.DEPEND_0.data))*1e-9;
    switch DATAC.tmMode
      case MMS_CONST.TmMode.srvy, error('kaboom')
      case  MMS_CONST.TmMode.slow, dtNominal = 20;
      case  MMS_CONST.TmMode.fast, dtNominal = 5;
      case  MMS_CONST.TmMode.brst, dtNominal = [0.625, 0.229];
      otherwise
        errS = 'Unrecognized tmMode';
        irf.log('critical',errS), error(errS)
    end
    flagOK = false;
    for i=1:numel(dtNominal)
      if dtSampling > dtNominal(i)*.95 && dtSampling < dtNominal(i)*1.05
        dtSampling = dtNominal(i); flagOK = true; break
      end
    end
    if ~flagOK
      errS = ['bad sampling for ' varPrefix probe '_enable'];
      irf.log('critical',errS), error(errS)
    end
    enabled.time = flag.DEPEND_0.data;
    nData = numel(enabled.time);
    enabled.data = zeros(nData,numel(fields));
    enabled.data(:,1) = flag.data;
    for iF=2:numel(fields)
      probe = fields{iF};
      flag = getv(dataObj,[varPrefix probe '_enable']);
      if isempty(flag)
        errS = ['cannot get ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      elseif numel(flag.data) ~= nData
        errS = ['bad size for ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      end
      enabled.data(:,iF) = flag.data;
    end
    newT = DATAC.(param).time;
    % Default to zero - probe disabled
    res = zeros(numel(newT), numel(fields));
    if all(diff(enabled.data))==0,
      ii = newT>enabled.time(1)-dtSampling & newT<=enabled.time(end);
      for iF=1:numel(fields), 
        res(ii,iF) = enabled.data(1,iF); 
      end
    else
      % TODO: implements some smart logic.
      errS = 'MMS_SDC_SDP_DATAMANAGER enabling/disabling probes not yet implemented.';
      irf.log('critical', errS); error(errS);
    end
  end
end


% Short function for verifying Time is increasing.
function check_monoton_timeincrease(time, dataType)
    
if(any(diff(time)<=0))
        err_str = ['Time is NOT increasing for the datatype ', dataType];
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:TIME:NONMONOTON', err_str);
end

end


% Short function for verifying variable is not stuck. DO NOT USE YET AS
% PRELIMINARY CDF FILES ARE STUCK.
function check_stuck_variable(var, varName)

% FIXME: Some values should perhaps be allowed to be the same for a limited
% number of datapoints, but for now check if ALL points in time are fixed.
if( all(diff(var))==0 )
    err_str = ['MMS_SDC_SDP_DATAMANAGER Variable ' varName,...
        ' appears to be stuck.'];
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:VARIABLE:STUCK', err_str);
end

end