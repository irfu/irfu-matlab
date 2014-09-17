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
        irf.log('warning',['First argument was not a dataobj but a file,'...
            ' trying to load with dataobj that file: ', dataObj, ...
            ', and store its data as: ',param,'.']);
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
        sig = {'e12','e34','e56'};
        init_param(sig)
        for iSig=1:length(sig)
          if ~isProbeEnabled(sig{iSig})
            DATAC.(param).(sig{iSig}).data = DATAC.(param).(sig{iSig}).data*NaN;
          end
        end
        
      case('dcv')
        sig = {'v1','v2','v3','v4','v5','v6'};
        init_param(sig)
        
        p1_off = ~isProbeEnabled('v1');
        p2_off = ~isProbeEnabled('v2');
        p3_off = ~isProbeEnabled('v3');
        p4_off = ~isProbeEnabled('v4');
        p5_off = ~isProbeEnabled('v5');
        p6_off = ~isProbeEnabled('v6');
        
        if p1_off && p2_off
          %DATAC.(param).v1.data = DATAC.(param).v1.data*NaN;
          %DATAC.(param).v2.data = DATAC.(param).v1.data;
        elseif p1_off
          if isfield(DATAC.dce,'e12') && ~isempty(DATAC.dce.e12)
            % Compute 
            % TODO:  implement real computation instead of this
          else
            %DATAC.(param).v1.data = DATAC.(param).v1.data*NaN;
          end
        elseif p2_off
          if isfield(DATAC.dce,'e12') && ~isempty(DATAC.dce.e12)
            % Compute 
            % TODO implement real computation instead of this
          else
            %DATAC.(param).v1.data = DATAC.(param).v1.data*NaN;
          end
        end
        
        % TODO: implement similar for p3-6
        
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

  function init_param(fields)
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
    if isempty(fields), return, end
    for iField=1:length(fields)
      DATAC.(param).(fields{iField}) = struct(...
        'data',sensorData(:,iField), ...
        'bitmask',zeros(size(sensorData(:,iField))));
    end
  end

  function res = isProbeEnabled(probe)
    %flag = get_variable(dataObj,[varPrefix probe '_enable']);
    
    flag = dataObj.data.([varPrefix probe '_enable']).data;
    if ~all(diff(flag))==0
      err_str = 'MMS_SDC_SDP_DATAMANAGER enabling/disabling probes not yet implemented.';
      irf.log('critical', err_str);
      error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    res = flag(1);
  end
end

function res = resample_bitmask(time,bmTime,bmData)

end


% Short function for verifying Time is increasing.
function check_monoton_timeincrease(time, dataType)
    
if(any(diff(time)<=0))
        err_str = ['MMS_SDC_SDP_DATAMANAGER Time is NOT increasing for the datatype ', dataType];
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