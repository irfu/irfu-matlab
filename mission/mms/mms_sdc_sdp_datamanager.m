function [ VarOut ] = mms_sdc_sdp_datamanager( varargin )
% MMS_SDC_SDP_DATAMANAGER load data into variables and then retrive them as
% needed. FIXME: ADD MORE HELP. Includes two small verification functions
% to ensure time is increasing and variables are not stuck at one single
% value.

narginchk(1,2); % One argument is simply retreive variable, Two arguments
% store "dataObj" as "dataType".

global DataInMemory; % Here we store read data.
global MMS_CONST; % Here we have Bitmask values for probes used.


if(nargin==2)
    
    % Make sure first argument is a dataobject (of a read cdf file).
    if( isa(varargin{1},'dataobj') )
        dataObj= varargin{1};
        
    elseif( ischar(varargin{1}) && exists(varargin{1}, 'file') && ...
            ischar(varargin{2}) )
        % If it is not a read cdf file, is it an unread cdf file? Read it.
        irf.log('warning',['First argument was not a dataobj but a file,'...
            ' trying to load with dataobj that file: ', varargin{1}, ...
            ', and store its data as: ',varargin{2},'.']);
        dataObj = dataobj(varargin{1}, 'tint', 0, 'true');
    
    else
        err_str = 'MMS_SDC_SDP_DATAMANAGER unknown input arguments.';
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
        
    end
    
    
    if( ischar(varargin{2}) )
        dataType = lower( varargin{2} );
        % Get first variable of CDF file, named mmsX_sdp_Epoch or similar.
        timeVar = getv( dataObj, dataObj.vars{1,1});
        check_monoton_timeincrease(timeVar.data, dataType);
        
        if( isfield(DataInMemory, dataType) )
            % Error, Warning or Notice for replacing the data variable?
            msg_str = ['MMS_SDC_SDP_DATAMANAGER replacing existing ', ...
                'variable with new data in ' dataType];
            irf.log('notice', msg_str);
            % warning('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', msg_str);;
        end
        
        switch(dataType)
            case('dce')
                DataInMemory.dce = [];
                DataInMemory.dce.time = timeVar.data;
                % FIXME: ADD CHECKS TO SEE WHICH PROBES ARE IN USE.
                % Also change variable number to match actual variable
                % name! Is now, 'mmsX_sdp_dce_sensors'.
                tmp = getv(dataObj, dataObj.vars{5,1});
                % FIXME: apply check if stuck.
                %check_stuck_variable(tmp.data(:,1), 'dce.p12');
                %check_stuck_variable(tmp.data(:,2), 'dce.p34');
                DataInMemory.dce.p12 = tmp.data(:,1);
                DataInMemory.dce.p34 = tmp.data(:,2);
                DataInMemory.dce.p56 = tmp.data(:,3);

            case('dcv')
                DataInMemory.dcv = [];
                DataInMemory.dcv.time = timeVar.data;
                % FIXME: ADD CHECKS TO SEE WHICH PROBES ARE IN USE.
                % Also change variable number to match actual variable
                % name! Is now, 'mmsX_sdp_dce_sensors'.
                tmp = getv(dataObj, dataObj.vars{5,1});
                % FIXME: apply check if stuck.
                %check_stuck_variable(tmp.data(:,1), 'dcv.p1');
                %check_stuck_variable(tmp.data(:,2), 'dcv.p3');
                DataInMemory.dcv.p1 = tmp.data(:,1);
                DataInMemory.dcv.p3 = tmp.data(:,2);
                DataInMemory.dcv.p5 = tmp.data(:,3);

            case('sunpulse')
                DataInMemory.sunpulse = [];
                % Also change variable number to match actual variable
                % name! Is now, 'mmsX_101_sunpulse'.
                tmp = getv( dataObj, dataObj.vars{19,1} );
                DataInMemory.sunpulse.timestamp = tmp.data;
                
            otherwise
                % Not yet implemented.
                err_str = ['MMS_SDC_SDP_DATAMANAGER unknown second ', ...
                    'parameter. The datatype ', dataType, ...
                    ' is not yet implemented'];
                irf.log('critical',err_str);
                error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
        
        end
        
    else
        err_str = 'MMS_SDC_SDP_DATAMANAGER second argument was not a string.';
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    
    end
    
end


if(nargin==1)
    % Should be name of Variable to return.
    if( ischar(varargin{1}) )
        varName = lower( varargin{1} );
        
        switch( varName )
            case('phase')
                % Phase, from sunpulse for now.
                if( isfield(DataInMemory.sunpulse,'phase') && ...
                        ~isempty(DataInMemory.sunpulse.phase) )
                    VarOut = DataInMemory.sunpulse.phase;
                else
                    % Calculate it, store it and return variable.
                    DataInMemory.sunpulse.phase = mms_sdc_sdp_phase( ...
                        DataInMemory.dce.time, ...
                        DataInMemory.sunpulse.timestamp);
                    VarOut = DataInMemory.sunpulse.phase;
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
                    err_str = ['The requested variable ', varName, ...
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
                    err_str = ['The requested variable ', varName, ...
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
    else
        err_str = 'MMS_SDC_SDP_DATAMANAGER unknown input argument. If only one argument is given it should be a string.';
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    
end

end


% Short function for verifying Time is increasing.
function check_monoton_timeincrease(time, dataType)
    
if(any(diff(time)<=0))
        err_str = ['MMS_SDC_SDP_DATAMANAGER Time is not monotonically increasing for the datatype ', dataType];
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:TIME:NONMONOTON', err_str);
end

end


% Short function for verifying variable is not stuck. DO NOT USE YET AS
% PRELIMINARY CDF FILES ARE STUCK.
function check_stuck_variable(var, varName)

% FIXME: A NUMBER OF DIFF(VAR) must be stuck, not just a sinlge value. But
% for now source file have all(diff(var)) == 0.
if( all(diff(var))==0 )
    err_str = ['MMS_SDC_SDP_DATAMANAGER Variable ' varName,...
        ' appears to be stuck.'];
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:VARIABLE:STUCK', err_str);
end

end