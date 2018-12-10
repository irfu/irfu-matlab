% Settings - Singleton class for global settings/constants used by BICAS.
% 
% Class for storing 
% 1) settings (variable values) which could reasonably be set via some user interface (CLI, configuration file), and
% 2) settings and constants which could reasonably be printed for the user.
%
%
% CONCEPT
% =======
% Data/settings are stored as key-value pairs. Keys are strings and values can be strings or numbers.
% --
% A settings object progress through three phases, in order, and stays definitively in the last phase:
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
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
classdef settings < handle
% BOGIQ: 
% ------
% PROPOSAL: Rename SIMPLE_DEMUXER. SIMPLE_CALIBRATION?
% PROPOSAL: Add extra information for every setting (key-value pair).
%   Ex: Human-readable description!
%   Ex: MATLAB class (data type)
%   Ex: Default value (so can display it if overridden)
%   Ex: Flag for origin of current value: default, config file, CLI argument.
%   Ex: Flag for write-protection (always use default value).
%       NOTE: Some settings (one?) make no sense to modify: config file path, STDOUT_PREFIX.
%   Ex: Flag for values which have not been set but must later be set.
%       Ex: MATLAB_COMMAND
%           CON: Is not really needed by BICAS.
%   --
%   NOTE: This information should only be given once in the code, and be hard-coded.
%
% PROPOSAL: Convention for "empty"/"not set"?!
%   TODO-DECISION/CON: Not really needed? Depends too much on the variable/setting.
%
% PROPOSAL: Move out interpretation of strings as numeric values?!
    
    properties(Access=private)
        defineDisabledForever = false;   % Whether defining new keys is disallowed or not. True if readOnlyForever==true.
        readOnlyForever       = false;   % Whether modifying the object is allowed or not.
        map;                             % Map containing settings data.
    end
    
    %###################################################################################################################
    
    methods(Access=public)
        
        % Constructor
        function obj = settings()
            % IMPLEMENTATION NOTE: "map" reset here since empirically it is not reset every time an instance is created
            % if it is only reset in the "properties" section. Otherwise the value from the previous execution is used
            % for unknown reasons.
            obj.map = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end
        
        
        
        function disable_define(obj)
            obj.defineDisabledForever = true;
        end
        
        
        
        function make_read_only(obj)
            obj.disable_define();
            obj.readOnlyForever = true;
        end
        
        
        
        % Define a NEW key and set the value.
        function define_setting(obj, key, value)
            % ASSERTIONS
            if obj.defineDisabledForever
                error('BICAS:settings:Assertion', 'Trying to define new keys in settings object which disallows defining new keys.')
            end
            if obj.map.isKey(key)
                error('BICAS:settings:Assertion', 'Trying to define pre-existing settings key.')
            end
            
            obj.map(key) = value;
        end
        
        
        
        % Set a pre-existing key value.
        function set_prexisting(obj, key, newValue)
            % ASSERTIONS
            if obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Trying to modify read-only settings object.')
            end
            if ~obj.map.isKey(key)
                error('BICAS:settings:Assertion', 'Trying to define non-existing settings key.')
            end
            
            oldValue = obj.map(key);
            
            if isnumeric(oldValue) && isnumeric(newValue)
                obj.map(key) = newValue;
            elseif ischar(oldValue) && ischar(newValue)
                obj.map(key) = newValue;
            else
                error('BICAS:settings:Assertion:IllegalArgument', ...
                    'New settings value either (1) does not match the type of the old settings value, or (2) is neither numeric nor char.')
            end

            obj.map(key) = newValue;
        end



        % Modify multiple settings, where the values are strings but converted to numerics as needed. Primarily intended
        % for updating settings with values from CLI arguments (which by their nature are initially strings).
        %
        %
        % ARGUMENTS
        % =========
        % ModifiedSettingsAsStrings : containers.Map
        %   <keys>   = Settings keys (strings). Must pre-exist as a SETTINGS key.
        %   <values> = Settings values AS STRINGS.
        %              Preserves the type of settings value for strings and numerics. If the pre-existing value is
        %              numeric, then the argument value will be converted to a number.
        %              Numeric row vectors are represented as a comma separated-list (no brackets), e.g. "1,2,3".
        %              Empty numeric vectors can not be represented.
        %
        % 
        % NOTE/BUG: No good checking (assertion) of whether the string format of a vector makes sense.
        % NOTE: Does not check if numeric vectors have the same size as old value.
        %
        function set_preexisting_from_strings(obj, ModifiedSettingsMap)
            
            % ASSERTION
            if obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Not allowed to modify read-only settings object.')
            end
            
            keysList = ModifiedSettingsMap.keys;
            for iModifSetting = 1:numel(keysList)
                key              = keysList{iModifSetting};
                newValueAsString = ModifiedSettingsMap(key);
                
                % ASSERTION
                if ~isa(newValueAsString, 'char')
                    error('BICAS:settings:Assertion:IllegalArgument', 'Map value is not a string.')
                end
                
                % Use old value to convert string value to appropriate MATLAB class.
                oldValue = obj.get_tv(key);   % ASSERTS: Key pre-exists
                if isnumeric(oldValue)
                    %newValue = str2double(newValueAsString);
                    newValue = textscan(newValueAsString, '%f', 'Delimiter', ',');
                    newValue = newValue{1}';   % Row vector.
                elseif ischar(oldValue)
                    newValue = newValueAsString;
                else
                    error('BICAS:settings:Assertion:ConfigurationBug', 'Can not handle the MATLAB class=%s of internal setting "%s".', class(oldValue), key)
                end
            
                % Overwrite old setting.
                obj.set_prexisting(key, newValue);
            end

        end



        % Return settings value for a given, existing key.
        % Only works when object is read-only, and the settings have their final values.
        %
        % IMPLEMENTATION NOTE: Short function name since function is called many times, often repeatedly.
        % FV = Final value
        function value = get_fv(obj, key)
            % ASSERTIONS
            if ~obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Not allowed to call this method for non-read-only settings object.')
            end
            if ~obj.map.isKey(key)
                error('BICAS:settings:Assertion:IllegalArgument', 'There is no setting "%s".', key)
            end
            
            value = obj.map(key);
        end
        
        
        
        % Return settings value for a given, existing key.
        % Only works when object is NOT read-only, and the settings might change their value later.
        %
        % IMPLEMENTATION NOTE: Name analogous to get_fv. Hence it is short.
        % TV = Tentative value
        function value = get_tv(obj, key)
            % ASSERTIONS
            if obj.readOnlyForever
                error('BICAS:settings:Assertion', 'Not allowed to call this method for read-only settings object.')
            end
            if ~obj.map.isKey(key)
                error('BICAS:settings:Assertion:IllegalArgument', 'There is no setting "%s".', key)
            end            
            
            value = obj.map(key);
        end
        
        
        
        function keyList = get_keys(obj)
            keyList = obj.map.keys;
        end
        
    end    % methods(Access=public)
    
end

