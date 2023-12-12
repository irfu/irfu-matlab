%
% Singleton class for global settings/constants used by BICAS, and which could
% reasonably be set via some user interface (default values, configuration file,
% CLI).
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
% A settings object progresses through three phases, in order, and stays
% write-protected in the last phase:
% (1) After creation: New keys can be defined and set to their initial values.
% (2) Definition disabled: Can set/modify the values of pre-existing keys
% (3) Read-only: Can not modify the object at all. Can only read key values.
%     (Object can not leave this phase.)
% Separate get methods are used for phases (1)-(2), and (3) respectively to
% ensure that the caller knows whether it reads a final value or not.
% --
% RATIONALE: This concept makes it natural to, when possible, clearly and
% conclusively separate the writing (setting) and reading of settings. Ideally,
% we would want all the writing to be followed by all the reading, but in
% practice they overlap and there does not seem to be a way of avoiding it in
% BICAS. For those cases it is useful to be forced to use a different get method
% to highlight that the read value is tentative (which it may be during
% initialization).
% 
%
% NOTE
% ====
% The class stores all overriden values, not just the latest ones. This has not
% been taken advantage of yet, but is intended for better logging the sources of
% settings and how they override each other. /2020-01-23
%
%
% ~BUG POTENTIAL: Support for 1D cell arrays may not be completely implemented.
%   ~BUG: Does not currently support setting 0x0 vectors (requires e.g. 0x1).
%         Inconvenient when working with values from CLI arguments and log files
%         since less convenient to write a 0x1 or 1x0 literal?!
%       Ex: Have to write zeros(0,1) instead of []?
% 
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-22
%
classdef Settings < handle
% BOGIQ: 
% ------
% PROPOSAL: Add extra information for every setting (key-value pair).
%   PROPOSAL: Human-readable description!
%   PROPOSAL: MATLAB class (data type)
%   PROPOSAL: Default value (so can display it if overridden)
%   PROPOSAL: Flag for write-protection (always use default value).
%       NOTE: Some settings (one?) make no sense to modify: config file path.
%   PROPOSAL: Flag for values which have not been set but must later be set.
%       PROPOSAL: MATLAB_COMMAND
%           CON: Is not really needed by BICAS.
%   PROPOSAL: Legal alternatives.
%       PRO: Rapid feedback when using bad value. Does not require triggering
%            code.
%       PRO: Clear in code (bicas.create_default_BSO()).
%       CON: Might not be consistent with how the settings values are actually used in the code.
%            Duplicates that decision.
%       PROPOSAL: Submit function (value-->boolean) that specifies what is legal
%                 and not. Can have set of pre-defined functions.
%           TODO-NI: How relates to how values are converted to display strings?
%           TODO-NI: How relates to how values are converted from strings (config file, CLI argument)?
%           
%       PROPOSAL: String constants.
%       PROPOSAL: Value type (MATLAB class)
%           Ex: Logical
%           CON: Not necessary since initial/default value specifies it.
%   --
%   NOTE: This information should only be given once in the code, and be hard-coded.
%
% PROPOSAL: Convention for "empty"/"not set"?!
%   TODO-DEC/CON: Not really needed? Depends too much on the variable/setting.
%
% PROPOSAL: Initialize by submitting map.
%   PRO: Can remove methods define_setting, disable_define.
%   CON: Can not easily add metadata for every variable (in the future),
%        e.g. permitted values (data type/class, range).
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
%   PROBLEM: bicas.interpret_config_file() and bicas.interpret_CLI_args() must
%            then be able to return info on a setting being set multiple times,
%            and return that information. As of now (2021-08-19) they only
%            return the final setting.
%       PROPOSAL: Submit BSO to those functions.
%           CON: Automatic testing becomes harder. Harder to test returned value. Harder to submit varied BSO.
%       PROPOSAL: Return KVPL.
%           NOTE: KVPL only permits string values(?).
%
% PROPOSAL: Make it possible to load multiple config files. Subsequent log files override each other.
%   TODO-DEC: Should the internal order of --set and --config arguments matter? Should a --config override a previous
%                  --set?
%
% PROPOSAL: Enforce keeping same MC (instead of the locally define "value
%           type").
%   CON: Could not specify values inside of cell arrays.
%       NOTE: Cell arrays are noy yet really supported /2023-12-12).



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(Access=private)
        % Whether defining new keys is disallowed or not. Always true if
        % readOnlyForever==true.
        defineDisabledForever = false;   
        
        % Whether modifying the object is allowed or not.
        readOnlyForever       = false;
        
        % Map containing the actual settings data.
        SkvMap;
    end
    
    
    
    properties(Constant)
        VALUE_SOURCE_DEFAULT = 'default';
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % Constructor
        function obj = Settings()
            % IMPLEMENTATION NOTE: "SkvMap" reset here since empirically it is
            % not reset every time an instance is created if it is only reset in
            % the "properties" section. Otherwise the value from the previous
            % execution is used for unknown reasons.
            obj.SkvMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end



        function disable_define(obj)
            obj.defineDisabledForever = true;
        end



        function make_read_only(obj)
            obj.disable_define();
            obj.readOnlyForever = true;
        end



        % Define a NEW key and set the corresponding value.
        %
        % NOTE: Key values in the form of MATLAB values, i.e. NOT string that
        % need to be parsed. cf .override_values_from_strings().
        %
        function define_setting(obj, key, defaultValue)
            % ASSERTIONS
            if obj.defineDisabledForever
                error('BICAS:Assertion', ...
                    ['Trying to define new keys in settings object which', ...
                    ' disallows defining new keys.'])
            end
            if obj.SkvMap.isKey(key)
                error('BICAS:Assertion:ConfigurationBug', ...
                    'Trying to define pre-existing settings key.')
            end
            bicas.Settings.assert_legal_value(defaultValue)
            
            
            
            % NOTE: Needs to be able to handle cell-valued values.
            Skv = bicas.SettingsKeyValue(...
                defaultValue, ...
                bicas.Settings.VALUE_SOURCE_DEFAULT);
            obj.SkvMap(key) = Skv;
        end



        % Set a PRE-EXISTING key value (i.e. override the default at the very
        % least) using MATLAB values.
        %
        % NOTE: Does not check if numeric vectors have the same size as old
        % value.
        % IMPLEMENTATION NOTE: BICAS does not need this method to be public, but
        % it is useful for other code (manual test code) to be able to override
        % settings using MATLAB values.
        function override_value(obj, key, newValue, valueSource)
            
            % ASSERTIONS
            irf.assert.castring(valueSource)
            if obj.readOnlyForever
                error('BICAS:Assertion', ...
                    'Trying to modify read-only settings object.')
            end
            
            Skv = obj.get_SKV_private(key);
            
            % ASSERTION: Old and new value have the same value type.
            if ~strcmp(...
                    bicas.Settings.get_value_type(newValue), ...
                    obj.get_setting_value_type(key))
                
                error('BICAS:Assertion:IllegalArgument', ...
                    ['New settings value does not match the type of the', ...
                    ' old settings value for key "%s".'], ...
                    key)
            end

            obj.SkvMap(key) = Skv.override(newValue, valueSource);
        end

        
        
        % Override multiple settings, where the values are strings but converted
        % to numerical values as needed. Primarily intended for updating
        % settings with values from CLI arguments and/or config file (which by
        % their nature have string values).
        %
        % NOTE: Method is essentially a for loop around .override_value().
        %
        % NOTE: Indirectly specifies the syntax for string values which
        % represent non-string-valued settings.
        %
        % NOTE/BUG: No good checking (assertion) of whether the string format of
        % a vector makes sense.
        %
        %
        % ARGUMENTS
        % =========
        % ModifiedSettingsMap
        %       containers.Map
        %       <keys>   = Settings keys (strings). Must pre-exist as a SETTINGS
        %                  key.
        %       <values> = Settings values AS STRINGS.
        %                  Preserves the type of settings value for strings and
        %                  numerics. For example, if the pre-existing value is
        %                  numeric, then the argument value will be converted to
        %                  a number. Numeric row vectors are represented as a
        %                  comma separated-list (no brackets), e.g. "1,2,3".
        %                  Empty numeric vectors can not be represented.
        %
        function obj = override_values_from_strings(...
                obj, ModifiedSettingsMap, valueSource)

            keysList = ModifiedSettingsMap.keys;
            for iModifSetting = 1:numel(keysList)
                key              = keysList{iModifSetting};
                newValueAsString = ModifiedSettingsMap(key);

                %==================================================
                % Convert string value to appropriate MATLAB class.
                %==================================================
                value_type = obj.get_setting_value_type(key);
                newValue   = bicas.Settings.convert_str_to_value(...
                    value_type, newValueAsString);

                % Overwrite old setting.
                obj.override_value(key, newValue, valueSource);
            end
        end



        function keyList = get_keys(obj)
            keyList = obj.SkvMap.keys;
        end
        
        
        
        % Return the settings value (that is actually going to be used) for a
        % given, existing key. Only works when object is read-only, and the
        % settings have their final values.
        %
        % IMPLEMENTATION NOTE: Short function name since function is called many
        % times, often repeatedly.
        % FV = Final value
        %
        % RETURN VALUES
        % ==============
        % value
        %       The value of the setting.
        % key
        %       The name of the settings key, i.e. identical to the argument
        %       "key".
        %       IMPLEMENTATION NOTE: This is useful in code that tries to avoid
        %       hardcoding the key string too many times. That way, the key is
        %       hard-coded once (in the call to this method), and then
        %       simultaneously assigned to a variable that is then used in the
        %       vicinity for error/warning/log messages etc. It is the second
        %       return value so that it can be ignored when the caller does not
        %       need it.
        function [value, key] = get_fv(obj, key)
            % ASSERTIONS
            if ~obj.readOnlyForever
                error('BICAS:Assertion', ...
                    ['Not allowed to call this method for non-read-only', ...
                    ' settings object.'])
            end

            Skv   = obj.get_SKV(key);
            value = Skv.valuesCa{end};
        end
        
        
        
        % Return settings value for a given, existing key. Only works when
        % object is read-only, and the settings have their final values.
        %
        % IMPLEMENTATION NOTE: Short function name since function is called many
        % times, often repeatedly.
        function Skv = get_SKV(obj, key)
            % ASSERTIONS
            if ~obj.readOnlyForever
                error('BICAS:Assertion', ...
                    ['Not allowed to call this method for a non-read-only', ...
                    ' settings object.'])
            end
            if ~obj.SkvMap.isKey(key)
                error('BICAS:Assertion:IllegalArgument', ...
                    'There is no setting "%s".', key)
            end

            Skv = obj.SkvMap(key);
        end
        


    end    % methods(Access=public)
    
    
    
    %##########################
    %##########################
    % PRIVATE INSTANCE METHODS
    %##########################
    %##########################
    methods(Access=private)
        
        
        
        % Needs to be public so that caller can determine how to parse string,
        % e.g. parse to number.
        function valueType = get_setting_value_type(obj, key)
            Skv = obj.get_SKV_private(key);
            
            % NOTE: Always use default/first value.
            valueType = bicas.Settings.get_value_type(Skv.valuesCa{1});
        end



        % Return SKV for a given, existing key.
        %
        % NOTE: Does not check obj.readOnlyForever since method is used
        % internally for non-readonly objects. Can therefore not be replaced by
        % public method get_SKV().
        %
        % RATIONALE: Exists to give better error message when using an illegal
        % key, than just calling obj.SkvMap directly.
        function Skv = get_SKV_private(obj, key)
            % ASSERTIONS
            if ~obj.SkvMap.isKey(key)
                error('BICAS:Assertion:IllegalArgument', ...
                    'There is no setting "%s".', key)
            end
            
            Skv = obj.SkvMap(key);
        end
        
        
        
    end    % methods(Access=private)
    
    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Access=private, Static)
        
        
        
        % Defines what is a legal value.
        function assert_legal_value(value)
            if ischar(value)
                
                % Do nothing
                irf.assert.castring(value)
                
            elseif isnumeric(value) ...
                    || iscell(value) ...
                    || islogical(value)
                
                irf.assert.vector(value)
                
            else
                
                error('BICAS:Assertion:IllegalArgument', ...
                    'Argument "value" is illegal.')
            end
        end
        
        
        
        function value = convert_str_to_value(settingValueType, valueAsString)
            % ASSERTION
            if ~isa(valueAsString, 'char')
                error('BICAS:Assertion:IllegalArgument', ...
                    'Map value is not a string.')
            end

            switch(settingValueType)

                case 'numeric'

                    % NORMALIZE: Remove surrounding brackets [] if present
                    % ----------------------------------------------------
                    % NOTE: textscan() can not handle [].
                    % IMPLEMENTATION NOTE: In principle, one could state that
                    % the format for specifying an array (as a string) is
                    % "1,2,3", but in practice, it is natural to use "[1,2,3]"
                    % since that is the MATLAB syntax. An earlier implementation
                    % falsely interpreted "[1,2,3]" as an empty array (i.e.
                    % without raising error), due to not considering this.
                    if regexp(valueAsString, '^\[.*\]$')
                        valueAsString = valueAsString(2:end-1);
                    end
                    
                    if isempty(valueAsString)
                        value = zeros(1,0);    % Zero-length row vector.
                    else
                        % NOTE: textscan() can not handle empty string (not even
                        %       1x0).
                        value = textscan(valueAsString, '%f', 'Delimiter', ',');

                        % IMPLEMENTATION NOTE: Not sure if/when this error message
                        % can be triggered.
                        if ~isscalar(value)
                            error('BICAS:Assertion:IllegalArgument', ...
                                'Can not parse string "%s" as array or scalar.', valueAsString)
                        end

                        value = value{1}';    % Row vector.
                    end

                case 'logical'
                    if strcmpi(valueAsString, 'true')
                        value = true;
                    elseif strcmpi(valueAsString, 'false')
                        value = false;
                    else
                        error('BICAS:Assertion:IllegalArgument', ...
                            'Can not parse supposed logical settings value "%s".', ...
                            valueAsString)
                    end

                case 'string'
                    value = valueAsString;

                otherwise
                    error('BICAS:Assertion:IllegalArgument', ...
                        ['Can not interpret argument settingValueType="%s"'], ...
                        valueAsString)
            end
            
            bicas.Settings.assert_legal_value(value)
        end
        
        
        
        % NOTE: Returns one of its own categories. Does not return MC.
        function valueType = get_value_type(value)
            if isnumeric(value)
                valueType = 'numeric';
            elseif islogical(value)
                valueType = 'logical';
            elseif ischar(value)
                valueType = 'string';
            else
                error('BICAS:ConfigurationBug', ...
                    'Settings value (old or new) has an illegal MATLAB class.')
            end
        end
        
        
        
    end    % methods(Access=private, Static)
    
end
