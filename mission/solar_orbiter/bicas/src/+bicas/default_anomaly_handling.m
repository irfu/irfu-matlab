%
% Given a setting for how to handle an anomaly that has already occurred.
%
% The concept is to be able to make various forms of error/anomaly handling
% shorter and more consistent, and make the use of anomaly handling settings
% more consistent. This code should be able to either handle or assist in the
% common cases:
% (1) Give a warning, but otherwise ignore.
% (2) Trigger error because of anomaly (like an assertion).
% (3) Trigger error because of illegal anomaly handling setting
%     (assertion on setting value).
% (4) Print warning message when outside code mitigate the anomaly/uses
%     mitigation/workaround.
%     NOTE: There could be multiple forms of mitigation/workarounds.
%
%
% DEFINITIONS
% ===========
% Anomaly : Something that should in theory never happen (and for which there
%           could in principle be an assertion error), but which can in practice
%           be handled in multiple depending on a setting.
%           Ex: Anomalies: Bad input dataset metadata which processing does not
%               need.
%           Ex: Handling: Error, warning, workaround/mitigation
%
%
% EXAMPLE 1: Handle ERROR, mitigation/workaround, illegal value
% =============================================================
% anomalyDescrMsg = 'Description of anomaly.';
% [settingValue, settingKey] = Bso.get_fv(...
%     'OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
% switch(settingValue)
%     case 'USE_WORKAROUND_1'
%         bicas.default_anomaly_handling(L, settingValue, settingKey, ...
%           'other', anomalyDescrMsg)
%         L.log('warning', 'Description of mitigation/workaround 1.')
%         % Code for mitigating/workaround 1.
%
%     case 'USE_WORKAROUND_2'
%         bicas.default_anomaly_handling(L, settingValue, settingKey, ...
%           'other', anomalyDescrMsg)
%         L.log('warning', 'Description of mitigation/workaround 2.')
%         % Code for mitigating/workaround 2.
%
%     otherwise
%         bicas.default_anomaly_handling(L, settingValue, settingKey, ...
%           'E+illegal', ...
%           anomalyDescrMsg, 'BICAS:SWMProcessing')
% end
%
%
% EXAMPLE 2: Handle ERROR, WARNING, illegal value
% ===============================================
% anomalyDescrMsg = 'Description of anomaly.';
% [settingValue, settingKey] = Bso.get_fv(...
%     'OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
% bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
%     anomalyDescrMsg, 'BICAS:SWMProcessing')
%
%
% ARGUMENTS
% =========
% casesHandled          : String constant. Which setting values are being
%                         handled by this particular call.
% anomalyDescriptionMsg : String.
% errorId               : String. Optional for casesHandled == 'other'
% --
% NOTE: Order of settingValue, settingKey.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-26
%
function default_anomaly_handling(L, settingValue, settingKey, ...
        casesHandled, anomalyDescriptionMsg, errorId)
    
    % PROPOSAL: Accept one struct "AnomalyHandlingInfo" for settingValue/Key, msg, errorId).
    %   PROPOSAL: Combine with Bso.get_fvs which returns struct with key+value.
    %   PRO: Shorter to call function repeatedly.
    %   PROPOSAL: Optional to use.
    %       NOTE: Implies that casesHandled should be earlier in argument list.
    %
    % PROBLEM: Want multi-row messages for log file (easier to read), but
    %          single-row messages for Exception.message (easier to grep) (?)
    %   PROPOSAL: Supply multi-row message. Always log it. Convert to single-row
    %             message (remove linefeed) for Exception.message.
    %   PROPOSAL: Use multi-row messages for both log file and
    %             Exception.message. User should use grep -Ax.
    %
    % PROBLEM: How handle mitigation using library functions?
    %   Library functions should not include application-specific warning/error
    %   messages, log messages, or setting values.
    %   PROBLEM: How handle library functions with multiple mitigations in sequence?
    %
    %
    %
    % PROPOSAL: Model for mitigation using with library functions and mutually
    %           exclusive forms of mitigation
    %   switch(settingValue)
    %       case 'MITIGATION_1'
    %           [y2, didMitigation1] = do_something(..., 'use mitigation 1');
    %           if didMitigation1
    %               bicas.default_anomaly_handling(...'other'...)
    %               L.log('Did mitigation 1.')
    %           end
    %       case 'MITIGATION_2'
    %           [y2, didMitigation2] = do_something(..., 'use mitigation 2');
    %           if didMitigation2
    %               bicas.default_anomaly_handling(...'other'...)
    %               L.log('Did mitigation 2.')
    %           end
    %       otherwise
    %           [y2, detectedAnomaly] = do_something(..., 'no mitigation');
    %           if detectedAnomaly
    %               bicas.default_anomaly_handling(...'E+W+illegal'...)
    %           end
    %   end
    %   PROPOSAL: Set casesHandled, mitigation description in respective case
    %             statement and then log mitigation description and call
    %             bicas.default_anomaly_handling once.
    %       CON: Still has to distinguish mitigations and non-mitigations
    %   PROPOSAL: Include logging mitigation message in default_anomaly_handling.
    %
    %
    %
    % PROPOSAL: Model for mitigation using with library functions and mutually
    %           exclusive forms of mitigation
    %   switch(settingValue)
    %       case 'MITIGATION_1'
    %           [y2, detectedAnomaly] = do_something(..., 'use mitigation 1');
    %           % detectedAnomaly also covers whether mitigation 1 was done.
    %           mitigationDescrMsg = 'Mitigation 1 description.';
    %           casesHandled       = 'other';
    %
    %       case 'MITIGATION_2'
    %           [y2, detectedAnomaly] = do_something(..., 'use mitigation 2');
    %           % detectedAnomaly also covers whether mitigation 2 was done.
    %           mitigationDescrMsg = 'Mitigation 2 description.';
    %           casesHandled       = 'other';
    %
    %       otherwise
    %           [y2, detectedAnomaly] = do_something(..., 'no mitigation');
    %           casesHandled       = 'E+W+illegal';
    %   end
    %   if detectedAnomaly
    %       bicas.default_anomaly_handling(..., casesHandled, ...)
    %       L.log(mitigationDescrMsg)
    %   end
    
    
            
    assert(ischar(settingKey), 'Argument settingKey is not a string.')
    
    setting1RowMsg = sprintf(...
        ['The behaviour when encountering this anomaly is', ...
        ' determined by setting %s = "%s".'], ...
        settingKey, settingValue);   % 1 row.
    setting2RowMsg = sprintf(...
        ['The behaviour when encountering this anomaly/error is', ...
        ' determined by setting\n    %s = "%s"\n'], ...
        settingKey, settingValue);   % 2 rows.
    
    
    PREFIX   = 'ANOMALY: ';
    N_INDENT = 4;
    ILLEGAL_SETTING_MSG = ...
        ['The setting value is illegal.', ...
        ' Can therefore not handle the error/anomaly.'];
    
    anomalyDescriptionMsg = irf.str.indent(...
        anomalyDescriptionMsg, numel(PREFIX));
    anomalyDescriptionMsg(1:numel(PREFIX)) = PREFIX;

    switch(casesHandled)
        case 'other'
            L.log('warning', anomalyDescriptionMsg)
            logi( 'warning', setting2RowMsg)
            return    % NOTE: RETURN
            
        case 'E+illegal'
            handleWarning = 0;
            
        case {'E+W+illegal', 'W+E+illegal'}
            handleWarning = 1;
            
        otherwise
            error('BICAS:Assertion', ...
                'Illegal argument caseshandled="%s"', casesHandled)
    end
    
    
    
    switch(settingValue)
        case 'ERROR'
            L.log('error', anomalyDescriptionMsg)
            logi( 'error', setting2RowMsg)
            error(errorId, '%s %s', anomalyDescriptionMsg, setting1RowMsg)
            
        case 'WARNING'
            if handleWarning
                L.log('warning', anomalyDescriptionMsg)
                logi( 'warning', ...
                    'Ignoring this anomaly (at least in this part of the code).')
                logi( 'warning', setting2RowMsg)
            else
                handle_illegal_settingValue()
            end

        otherwise
            handle_illegal_settingValue()
    end
    
    %===========================================================================
    function handle_illegal_settingValue()
        L.log('error', anomalyDescriptionMsg)
        logi( 'error', setting2RowMsg)
        logi( 'error', ILLEGAL_SETTING_MSG)
        error('BICAS:ConfigurationBug', ...
            '%s %s %s', ...
            anomalyDescriptionMsg, setting1RowMsg, ILLEGAL_SETTING_MSG)
        % NOTE: Function defines its own ID.
    end
    %===========================================================================
    function logi(logLevel, str)
        L.log(logLevel, irf.str.indent(str, N_INDENT))
    end
end
