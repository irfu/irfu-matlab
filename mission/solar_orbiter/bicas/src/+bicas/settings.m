%
% Singleton class for global settings/constants used by BICAS, and which could reasonably be set via some
% user interface (default values, configuration file, CLI).
%
%
% CONCEPT
% =======
% Data/settings are stored as a set of key-value pairs.
%   Keys : String
%   Value : One of below:
%       (1) strings
%       (2) numbers (1D vector)
%       (3) cell array of strings (1D vector)
% --
% A settings object progress through three phases, in order, and stays ROC_PIP_NAME/write-protected in the last phase:
% (1) From creation: New keys can be defined and set to their initial values.
% (2) Definition disabled: Can set the values of pre-existing keys
% (3) Read-only: Can not modify the object at all. Can only read key values. (Object can not leave this phase.)
% Separate get methods are used for phases (1)-(2), and (3) respectively.
% RATIONALE: This concept makes it natural to, when possible, clearly and conclusively separate the writing (setting)
% and reading of settings. Ideally, we would want all the writing to be followed by all the reading, but in practice
% they overlap and there does not seem to be a way of avoiding it in BICAS. For those cases it is useful to be forced to
% use a different get method to highlight that the read value is tentative (which it may be during initialization).
% 
%
% NOTE
% ====
% Class stores all overriden values, not just the latest ones. This has not been taken advantage of yet, but is
% intended for better logging the sources of settings and how they override each other. /2020-01-23
%
%
% ~BUG POTENTIAL: Support for 1D cell arrays may not be completely implemented.
%   ~BUG: Does not currently support setting 0x0 vectors (requires e.g. 0x1). Inconvenient.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
classdef settings < handle
% BOGIQ: 
% ------
% PROPOSAL: Add extra information for every setting (key-value pair).
%   PROPOSAL: Human-readable description!
%   PROPOSAL: MATLAB class (data type)
%   PROPOSAL: Default value (so can display it if overridden)
%   PROPOSAL: Flag for write-protection (always use default value).
%       NOTE: Some settings (one?) make no sense to modify: config file path, STDOUT_PREFIX.
%   PROPOSAL: Flag for values which have not been set but must later be set.
%       PROPOSAL: MATLAB_COMMAND
%           CON: Is not really needed by BICAS.
%   PROPOSAL: Legal alternatives.
%       PROPOSAL: String constants.
%       PROPOSAL: Value type (MATLAB class)
%           Ex: Logical
%           CON: Not necessary since initial/default value specifies it.
%   --
%   NOTE: This information should only be given once in the code, and be hard-coded.
%
% PROPOSAL: Convention for "empty"/"not set"?!
%   TODO-DECISION/CON: Not really needed? Depends too much on the variable/setting.
%
% PROPOSAL: Initialize by submitting map.
%   PRO: Can remove methods define_setting, disable_define.
%   CON: Can not easily add metadata for every variable (in the future), e.g. permitted values (data type/class, range).
%
% PROPOSAL: Be able to make some settings (default values) write-protected, not overridable.
%   CON: Of limited value.
%
% PROPOSAL: Store which settings were invoked (read) during a run.
%   PRO: Can summarize (and log) which settings are actually being used.
%   CON: Must distinguish between retrieving settings for actual use in algorithm, or for just logging.
%
% PROPOSAL: Enable BICAS to log where a key is set, and how many times. To follow how default value is overridden, and
%           how it is overriden twice in the same interface (in the config file or CLI arguments)
%   Ex: Config file specifies a new "default" value which is then overridden further below.
%   PROBLEM: interpret_config_file and interpret_CLI_args must then be aware of setting a key multiple times, and return
%   that information.
%       PROPOSAL: Submit SETTINGS to those functions.
%           CON: Automatic testing becomes harder. Harder to test returned value. Harder to submit varied SETTINGS.
%       PROPOSAL: Return KVPL.
%           NOTE: KVPL only permits string values(?).
%
% PROPOSAL: Make it possible to load multiple config files. Subsequent log files override each other.
%   TODO-DECISION: Should the internal order of --set and --config arguments matter? Should a --config override a previous
%                  --set?


    properties(Access=private)
        defineDisabledForever = false;   % Whether defining new keys is disallowed or not. Always true if readOnlyForever==true.
        readOnlyForever       = false;   % Whether modifying the object is allowed or not.
        DataMap;                         % Map containing the actual settings data.
    end



    %###################################################################################################################

        
    
    methods(Access=public)



        % Constructor
        function obj = settings()
            % IMPLEMENTATION NOTE: "DataMap" reset here since empirically it is not reset every time an instance is created
            % if it is only reset in the "properties" section. Otherwise the value from the previous execution is used
            % for unknown reasons.
            obj.DataMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end



        function disable_define(obj)
            obj.defineDisabledForever = true;
        end



        function make_read_only(obj)
            obj.disable_define();
            obj.readOnlyForever = true;
        end



        % Define a NEW key and set the corresponding value.
        function define_setting(obj, key, defaultValue)
            % ASSERTIONS
            if obj.defineDisabledForever
                error('BICAS:settings:Assertion', 'Trying to define new keys in settings object which disallows defining new keys.')
            end
            if obj.DataMap.isKey(key)
                error('BICAS:settings:Assertion:ConfigurationBug', 'Trying to define pre-existing settings key.')
            end
            
            % ASSERTIONS
            if ischar(defaultValue)
                % Do nothing
            elseif isnumeric(defaultValue) || iscell(defaultValue)
                EJ_library.assert.vector(defaultValue)
            else
                error('BICAS:settings:Assertion:IllegalArgument', 'Argument defaultValue is illegal.')
            end
            
            
            
            % NOTE: Needs to be able to handle cell-valued values.
            Setting = struct(...
                'value',       {defaultValue}, ...
                'valueSource', {'default'});
            assert(isscalar(Setting))
            obj.DataMap(key) = Setting;
        end



        % Set a PRE-EXISTING key value (i.e. override the default at the very least).
        % NOTE: Does not check if numeric vectors have the same size as old value.
        function override_value(obj, key, newValue, valueSource)
            % NOTE: Used to be public method. Can/should probably be rewritten or merged with set_preexisting_from_strings.
            
            % ASSERTIONS
            EJ_library.assert.castring(valueSource)
            if obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Trying to modify read-only settings object.')
            end
            
            valueArayStruct = obj.get_value_array_struct(key);
            
            % ASSERTION
            if ~strcmp(bicas.settings.get_value_type(newValue), obj.get_setting_value_type(key))
                error('BICAS:settings:Assertion:IllegalArgument', ...
                    'New settings value does not match the type of the old settings value.')
            end

            % IMPLEMENTATION NOTE: obj.DataMap(key).value = newValue;   % Not permitted by MATLAB.
            valueArayStruct(end+1).value       = newValue;
            valueArayStruct(end  ).valueSource = valueSource;
            obj.DataMap(key) = valueArayStruct;
        end

        
        
        function keyList = get_keys(obj)
            keyList = obj.DataMap.keys;
        end
        
        
        
        % Return the settings value (that is actually going to be used) for a given, existing key.
        % Only works when object is read-only, and the settings have their final values.
        %
        % IMPLEMENTATION NOTE: Short function name since function is called many times, often repeatedly.
        % FV = Final value
        %
        % RETURN VALUES
        % ==============
        % value : The value of the setting.
        % key   : The name of the settings key, i.e. identical to the argument "key".
        %         IMPLEMENTATION NOTE: This is useful in code that tries to avoid hardcoding the key string too many
        %         times. That way, the key is hardcoded once (in the call to this method), and then simultaneously
        %         assigned to a variable that is then used in the vicinity for error/warning/log messages etc. It is the
        %         second return value so that it can be ignored when the caller does not need it.
        function [value, key] = get_fv(obj, key)
            % ASSERTIONS
            if ~obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Not allowed to call this method for non-read-only settings object.')
            end
            valueStructArray = obj.get_value_array_struct(key);

            value = valueStructArray(end).value;
        end
        
        
        
        % Return settings value for a given, existing key.
        % Only works when object is read-only, and the settings have their final values.
        %
        % IMPLEMENTATION NOTE: Short function name since function is called many times, often repeatedly.
        % FV = Final value
        function valueArrayStruct = get_final_value_array(obj, key)
            % ASSERTIONS
            if ~obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Not allowed to call this method for non-read-only settings object.')
            end
            if ~obj.DataMap.isKey(key)
                error('BICAS:settings:Assertion:IllegalArgument', 'There is no setting "%s".', key)
            end
            
            
            valueArrayStruct = obj.DataMap(key);
            EJ_library.assert.struct(valueArrayStruct, {'value', 'valueSource'}, {})
        end
        


        % Needs to be public so that caller can determine how to parse string, e.g. parse to number.
        function valueType = get_setting_value_type(obj, key)
            valueArrayStruct = obj.get_value_array_struct(key);            
            valueType        = bicas.settings.get_value_type(valueArrayStruct(1).value);    % NOTE: Always use default/first value.
        end



    end    % methods(Access=public)
    
    
    
    methods(Access=private)
        
        
        
        % Return settings array struct for a given, existing key.
        %
        % RATIONALE: Exists to give better error message when using an illegal key, than just calling obj.DataMap
        % directly.
        function S = get_value_array_struct(obj, key)
            % ASSERTIONS
            if ~obj.DataMap.isKey(key)
                error('BICAS:settings:Assertion:IllegalArgument', 'There is no setting "%s".', key)
            end
            
            S = obj.DataMap(key);
        end
        
        
        
    end    % methods(Access=private)
    
    
    
    methods(Access=private, Static)
        
        
        
        function valueType = get_value_type(value)
            if isnumeric(value)
                valueType = 'numeric';
            elseif ischar(value)
                valueType = 'string';
            else
                error('BICAS:settings:ConfigurationBug', 'Settings value (old or new) has an illegal MATLAB class.')
            end
        end
        
        
        
    end
    
end
