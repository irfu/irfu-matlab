function [ VarOut ] = mms_sdc_sdp_datamanager( param, dataObj )
% mms_sdc_sdp_datamanager will store and retrive data for
%  	mms_sdc_sdp_datamanager( dataType, dataObj ) will store
%  	appropriate data variables related to dataType in the global variable
%  	DataInMemory for later retrival.
%  
%   [varOut] = mms_sdc_sdp_datamanager( variable ) will return the variable
%   requested to varOut, if no such variable has been created already it
%   will try and calculate it using the stored data.
%   
%  	Example:
%.scId
%  		mms_sdc_sdp_datamanager('dce',dceDataObj)
%  		mms_sdc_sdp_datamanager('phase')
%  
%   	See also DATAOBJ, MMS_SDC_SDP_CDF_IN_PROCESS.

narginchk(1,2); % One argument is simply retreive variable, Two arguments
% store "dataObj" as "dataType".

global DataInMemory; % Here we store read data.
global MMS_CONST; % Here we have Bitmask values for probes used.

if ~ischar(param),
  err_str = 'PARAM must be a string';
  irf.log('critical', err_str);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
end
if strcmpi(param, 'init')
  % Initialize
  if ~isnumeric(dataObj) || isempty(intersect(1:4,dataObj))
    err_str = 'second argument for INIT must be MMS_ID 1..4';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
  end
  DataInMemory = [];
  DataInMemory.scId = dataObj;
  DataInMemory.dce = [];
  DataInMemory.dcv = [];
  DataInMemory.hk_101 = [];
  return
end
if ~isfield(DataInMemory, 'scId')
  err_str = 'Data mamager not initialized! Run: mms_sdc_sdp_datamanager(''init'',scId)';
  irf.log('critical', err_str);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
end
param = lower(param);

if(nargin==2)   
    % Make sure first argument is a dataobject (of a read cdf file).
    if isa(dataObj,'dataobj') % do nothing
    elseif ischar(dataObj) && exists(dataObj, 'file')
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
    
    if( isfield(DataInMemory, param) ) && ~isempty(DataInMemory.(param))
      % Error, Warning or Notice for replacing the data variable?
      err_str = ['MMS_SDC_SDP_DATAMANAGER replacing existing ', ...
        'variable with new data in ' param];
      irf.log('critical', err_str);
      error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    
    varPrefix = sprintf('mms%d_sdp_',DataInMemory.scId);
    
    switch(param)
      case('dce')
        sig = {'e12','e34','e56'};
        init_param(sig)
        for iSig=1:length(sig)
          if ~isProbeEnabled(sig{iSig})
            DataInMemory.(param).(sig{iSig}).data = DataInMemory.(param).(sig{iSig}).data*NaN;
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
          %DataInMemory.(param).v1.data = DataInMemory.(param).v1.data*NaN;
          %DataInMemory.(param).v2.data = DataInMemory.(param).v1.data;
        elseif p1_off
          if isfield(DataInMemory.dce,'e12') && ~isempty(DataInMemory.dce.e12)
            % Compute 
            % TODO:  implement real computation instead of this
          else
            %DataInMemory.(param).v1.data = DataInMemory.(param).v1.data*NaN;
          end
        elseif p2_off
          if isfield(DataInMemory.dce,'e12') && ~isempty(DataInMemory.dce.e12)
            % Compute 
            % TODO implement real computation instead of this
          else
            %DataInMemory.(param).v1.data = DataInMemory.(param).v1.data*NaN;
          end
        end
        
        % TODO: implement similar for p3-6
        
      case('hk_101')
        varPrefix = sprintf('mms%d_101_',DataInMemory.scId);
        DataInMemory.(param) = [];
        DataInMemory.(param).dataObj = dataObj;
        x = getdep(dataObj,[varPrefix 'cmdexec']);
        DataInMemory.(param).time = x.DEPEND_O;
        check_monoton_timeincrease(DataInMemory.(param).time, param);
        % Add sunpulse times (TT2000) of last recieved sunpulse.
        DataInMemory.(param).sunpulse = dataObj.data.([varPrefix 'sunpulse']).data;
        % Add sunpulse indicator, real: 0, SC pseudo: 1, CIDP pseudo: 2.
        DataInMemory.(param).sunsps = dataObj.data.([varPrefix 'sunssps']).data;
        % Add CIDP sun period (in microseconds, 0 if sun pulse not real.
        DataInMemory.(param).iifsunper = dataObj.data([varPrefix 'iifsunper']).data;
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
                if isfield(DataInMemory,'dce') && ...
                    isfield(DataInMemory.dce,'phase') && ...
                    ~isempty(DataInMemory.dce.phase)
                  VarOut = DataInMemory.dce.phase;
                else
                    % Calculate it, store it and return variable.
                    DataInMemory.dce.phase = mms_sdc_sdp_phase( ...
                        DataInMemory.dce.time, ...
                        DataInMemory.hk_101.sunpulse);
                    VarOut = DataInMemory.dce.phase;
                end
                
                
                %%
                %THIS SHOULD BE CHANGED BUT FOR TESTING AND VERIFICATION
                %PURPOSES IT IS DONE THIS WAY. THE FOLLOWING MUST BE
                %CHANGED.
            case('dcetime')
                % Timestamp of dce
                if( isfield(DataInMemory.dce, 'time') && ...
                        ~isempty(DataInMemory.dce.time) )
                    VarOut = DataInMemory.dce.time;
                else
                    % Error
                    err_str = ['The requested variable ', param, ...
                        'does not exist.'];
                    irf.log('critical',err_str);
                    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
                end
                
            case('dcvtime')
                % Timestamp of dce
                if( isfield(DataInMemory.dcv, 'time') && ...
                        ~isempty(DataInMemory.dcv.time) )
                    VarOut = DataInMemory.dcv.time;
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
    DataInMemory.(param) = [];
    if ~all(diff(dataObj.data.([varPrefix 'samplerate_' param]).data)==0)
      err_str = 'MMS_SDC_SDP_DATAMANAGER changing sampling rate not yet implemented.';
      irf.log('warning', err_str);
      %error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    DataInMemory.(param).dataObj = dataObj;
    x = getdep(dataObj,[varPrefix param '_sensor']);
    DataInMemory.(param).time = x.DEPEND_O;
    check_monoton_timeincrease(DataInMemory.(param).time, param);
    sensorData = dataObj.data.([varPrefix param '_sensor']).data;
    if isempty(fields), return, end
    for iField=1:length(fields)
      DataInMemory.(param).(fields{iField}) = struct(...
        'data',sensorData(:,iField), ...
        'bitmask',zeros(size(sensorData(:,iField))));
    end
  end

  function res = isProbeEnabled(probe)
    flag = dataObj.data.([varPrefix probe '_enable']).data;
    if ~all(diff(flag))==0
      err_str = 'MMS_SDC_SDP_DATAMANAGER enabling/disabling probes not yet implemented.';
      irf.log('critical', err_str);
      error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    res = flag(1);
  end
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