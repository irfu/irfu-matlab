%
% Class which encapsulates the information stored for one settings key in
% bicas.Settings.
%
%
% IMPLEMENTATION NOTE
% ===================
% The class stores all overriden values, not just the latest one. This has not
% yet been taken advantage of, but is intended for making possible to do better
% logging of the sources of settings and how they override each other.
% /2023-12-13
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SettingsKeyValue
% PROPOSAL: Better name
%   ~Settings
%       CON: Want to avoid if moving class to package for settings.
%   ~Entry
%   ~Key-value pair
%   Value
%       CON: There are multiple stored values.
%           CON: Only one is really used.
%   Key
%       PRO: There is one key per object but multiple values.
%   SettingsEntry
%       CON: Two-letter abbreviation.
%   SettingsKeyEntry = SKE
%       CON: SKE is bad abbreviation. Foud in too many works, variables.
%   SettingsKeyValueEntry = SKVE
%       NOTE: Can drop "Settings" if moving class to settings package.
%
% PROPOSAL: Field for valueType.
% PROPOSAL: Enforce keeping same MC (instead of the locally defined "value
%           type").
%   CON: Could not specify values inside of cell arrays.
%       NOTE: Cell arrays are not yet really supported /2023-12-12).
%       CON: Can not by extending curernt scheme either.
%   CON: Requires more complex conversion string-->MATLAB value.
%       Ex: How specify type of integer as string?
%
% PROPOSAL: Add extra information for every setting (key-value pair).
%   PROPOSAL: Human-readable description!
%   PROPOSAL: Flag for values which have not been set but must later be set.
%       PROPOSAL: MATLAB_COMMAND
%           CON: Is not really needed by BICAS.
%   PROPOSAL: Legal values.
%       Ex: Legal string constants.
%       --
%       PRO: Rapid feedback when using bad value. Does not require triggering
%            code that actually uses value.
%       PRO: Clear in code (bicas.create_default_BSO()).
%       CON: Might not be consistent with how the settings values are actually
%            used in the code. Duplicates that decision.
%       PROPOSAL: Submit function handle (value-->boolean) that specifies what
%                 is legal and not. Can have set of pre-defined functions.
%           TODO-NI: How relates to how values are converted to display strings?
%           TODO-NI: How relates to how values are converted from strings (config file, CLI argument)?
%
%       PROPOSAL: Value type (MATLAB class)
%           Ex: Logical
%           CON: Not necessary since initial/default value specifies it.
%   --
%   NOTE: This information should only be given once in the code, and be hard-coded.
%
% PROPOSAL: Be able to make some settings (default values) write-protected, not overridable.
%   Ex: MATLAB_COMMAND
%   Ex: Config file path.
%   CON: Of limited value.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=public)
        % Column cell array of values.
        %
        % NOTE: Must be CA of values since the values might (1) be strings, or
        % (2) arrays of varying sizes.
        valuesCa

        % Column cell array.
        valueSourcesCa
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = SettingsKeyValue(initialValue, valueSource)
            [~] = bicas.SettingsKeyValue.assert_legal_value_get_type(initialValue);
            assert(ischar(valueSource))

            obj.valuesCa       = {initialValue};
            obj.valueSourcesCa = {valueSource};
        end



        function obj = override(obj, newValue, valueSource)
            assert(ischar(valueSource))
            assert(~ismember(valueSource, obj.valueSourcesCa))

            oldValueType = obj.get_value_type();
            newValueType = bicas.SettingsKeyValue.assert_legal_value_get_type(newValue);
            if ~strcmp(oldValueType, newValueType)
                error('BICAS:IllegalOverridingSettingValueType', ...
                    'Overriding setting will illegal value type.')
            end

            % NOTE: Creating column arrays.
            obj.valuesCa{end+1, 1}       = newValue;
            obj.valueSourcesCa{end+1, 1} = valueSource;
        end



        function n = N_values(obj)
            n = numel(obj.valuesCa);
        end



        % NOTE: Is a public function so that
        % bicas.Settings.override_values_from_strings() can use value type to
        % convert strings to MATLAB values.
        function valueType = get_value_type(obj)
            value     = obj.valuesCa{1};
            valueType = bicas.SettingsKeyValue.assert_legal_value_get_type(value);
        end



    end    % methods(Access=public)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % (1) Verifies that value is legal,
        % (2) returns type of value (so that it can be checked against new
        %     value if overriding it).
        %     NOTE: Returns one of its own categories. Does not return MC.
        function valueType = assert_legal_value_get_type(value)
            if ischar(value)
                irf.assert.castring(value)
                valueType = 'string';
            elseif isnumeric(value)
                irf.assert.vector(value)
                valueType = 'numeric';
            elseif islogical(value)
                irf.assert.vector(value)
                valueType = 'logical';
            else
                error('BICAS:Assertion:IllegalArgument', 'Argument "value" is illegal.')
            end
        end



    end



end
