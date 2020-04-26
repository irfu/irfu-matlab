%
% Given a setting for how to handle an anomaly that has already occurred, handle standardized setting values.
% Anomaly = something that should not happen but may, depending on a setting, be handled in multiple ways.
%
% The concept is to be able to make various forms of error/anomaly handling shorter and more consistent, and make the
% use of anomaly handling settings more consistent. This code should handle the common cases:
% (1) Give a warning, but otherwise ignore.
% (2) Give error because of anomaly.
% (3) Error because of illegal anomaly handling setting.
% (4) Mitigate the anomaly (print warning message)
%   The actual mitigation code is outside this function.
%  NOTE: The there could be multiple forms of mitigation.
%
%
% EXAMPLE USE
% ===========
% anomalyDescrMsg = '...';
% [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
% switch(settingValue)
%     case 'USE_FILLVAL'
%         nonmitigation_anomaly_handling(L, settingValue, settingKey, 'other', ...
%           anomalyDescrMsg)
%         L.log('warning', 'Mitigation description')
%         % Mitigate
%     otherwise
%         nonmitigation_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
%           anomalyDescrMsg, 'BICAS:execute_sw_mode:SWModeProcessing')
% end
%
%
% ARGUMENTS
% =========
% casesHandled          : String constant. Which setting values are being handled by this particular call.
% anomalyDescriptionMsg : String.
% errorId               : String. Optional for casesHandled == 'other'
% --
% NOTE: Order of settingValue, settingKey.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-26
%
function default_anomaly_handling(L, settingValue, settingKey, casesHandled, anomalyDescriptionMsg, errorId)

    assert(ischar(settingValue))
    assert(ischar(settingKey))
    
    %[settingValue, settingKey] = SETTINGS.get_fv(settingKey);
    
    setting1Msg = sprintf('The behaviour when encountering this error/anomaly is determined by setting %s = "%s".', ...
        settingKey, settingValue);   % 1 row.
    setting2Msg = sprintf('The behaviour when encountering this error/anomaly is determined by setting\n    %s = "%s"\n', ...
        settingKey, settingValue);   % 2 rows.
    
    

    switch(casesHandled)
        case 'other'
            L.log('warning', anomalyDescriptionMsg)
            L.log('warning', setting2Msg)
            return    % NOTE: RETURN
            
        case 'E+illegal'
            handleWarning = 0;
            
        case {'E+W+illegal', 'W+E+illegal'}
            handleWarning = 1;
            
        otherwise
            error('BICAS:default_anomaly_handling:Assertion', 'Illegal argument caseshandled="%s"', casesHandled)
    end
    
    
    
    switch(settingValue)
        case 'ERROR'
            L.log('error', anomalyDescriptionMsg)
            L.log('error', setting2Msg)
            error(errorId, '%s %s', anomalyDescriptionMsg, setting1Msg)
            
        case 'WARNING'
            if handleWarning
                L.log('warning', anomalyDescriptionMsg)
                L.log('warning', setting2Msg)
            else
                handle_illegal_settingValue()
            end

        otherwise
            handle_illegal_settingValue()
    end
    
    %=================================================================================================================
    function handle_illegal_settingValue()
        illegalSettingMsg = 'The setting value is illegal. Can therefore not handle the error/anomaly.';
        
        L.log('error', anomalyDescriptionMsg)
        L.log('error', setting2Msg)
        L.log('error', illegalSettingMsg)
        error('BICAS:Assertion:ConfigurationBug', ...
            '%s %s %s', anomalyDescriptionMsg, setting1Msg, illegalSettingMsg)
    end
end
