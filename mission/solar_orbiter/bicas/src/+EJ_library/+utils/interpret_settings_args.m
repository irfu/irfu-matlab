%
% Utility function for interpreting settings argument(s) on a standardized syntax/format passed to functions. 
%
%
% INTENDED USAGE
% ==============
% The function interprets a 1D cell array, representing a list of arguments to another function, typically the
% external function that calls this function, and typically the external function's last arguments by using varargin.
% --
% In order to require/permit a certain set of settings (keys), use an assertion for the set of field names in the
% returned struct. Note that if there are default values for all keys, then DefaultSettings can be used to obtain the
% list of field names.
%
%
% ALGORITHM
% =========
% The algorithm uses/constructs three analogous "Settings" structs.
% DefaultSettings      = argument.
% SettingsArg1         = argList{1}, if argList{1} exists and is a struct. Otherwise empty struct.
% SettingsArgListPairs = remainder of argList (excluding SettingsArg1), interpreted as pairs of field name + field value.
% --
% The functions returns a struct which is a combination (union of fields) of 
% (1) SettingsArgListPairs
% (2) SettingsArg1
% (3) DefaultSettings
% where field values are taken from the first top-most struct with that field, i.e. a higher one has precedence over a
% lower one, e.g. (1) has precedence over (2). Therefore, the structs do not need to, but are allowed to, contain the same fields.
%
%
% LIMITATIONS
% ===========
% Can not handle situations where set of fields in the "full settings struct" depends on other fields.
%   Ex: Settings.enabled = 1 implies there should be a field Settings.parameter, but Settings.enabled = 0 implies there
%       should not.
%
%
% ARGUMENTS
% =========
% DefaultSettings : Struct with default Settings (to be processed; not Settings for this functions).
% argList         : Cell array representing a sequence of arguments (presumably "varargin" or subset thereof) from
%                   another function that uses this function.
%                   It is either
%                       (1) {              key1,value1, ..., keyN,valueN}
%                       (2) {SettingsArg1, key1,value1, ..., keyN,valueN}
%                   NOTE: The cell array is permitted to be empty.
%
%
% RETURN VALUES
% =============
% Settings    : Struct. See algorithm.
% argList     : Cell array of strings, representing list of arguments passed to other function. Typically varargin
%               as received from the enclosing function directly.
%
%
% Initially created 2018-07-18 by Erik P G Johansson.
%
function [Settings] = interpret_settings_args(DefaultSettings, argList)
    % PROPOSAL: Assert SettingsArg1 and SettingsArgListPairs fields to always exist in DefaultSettings.
    %   CON: Makes it impossible to have default values for some settings/fields, but not for others.
    %   PROPOSAL: Option/flag for this behaviour.
    %   NOTE: User can easily add an assertion after:
    %       EJ_library.utils.assert.struct2(Settings, fieldnames(DEFAULT_SETTINGS), {})
    %
    % PROPOSAL: Reorg algorithm to produce better error messages when using bad combinations of arguments.
    %
    % PROPOSAL: Argument for required settings fields (DefaultSettings contains the ones which are optional in varargin).
    %   CON: Might have situations where the required settings fields depend on the the value of other settings fields.
    %   PROPOSAL: User should call EJ_library.utils.assert.struct2 instead. This is almost as succint.
    %
    % PROPOSAL: Add ability to recognize "string keywords" (one argument, instead of keyword+value) which indicate that
    % a flag (false/true) shall be set.
    %   NOTE: The default value is always false.
    %   TODO-DECISION: How represent in a struct that a field represents an optional string keyword?
    %       PROPOSAL: Special value, e.g. "string keyword", "argument keyword".
    %       PROPOSAL: Field name naming convention.
    %       NOTE: Might want it to be possible to both specify either a string keyword or a setting+value for the same
    %               setting e.g. if passing on variable values.
    %       PROBLEM: Does not want to restrict string keywords to possible field names (e.g. no whitespace)?
    %   CON: Can just as well just use key+value, and set value to 0/1. It is quite short.
    %
    % PROPOSAL: Use containers.Map instead of struct.
    %   PRO: Can have arbitrary keys with whitespace and non-text characters.
    %       Ex: "Extrapolate Z(omega>0) to Z(0)')"
    %   CON: Has less functionality for working with containers.Map than for struct. Union, difference, overwrite overlap etc.
    %       Ex: Want to internally merge/overwrite different sources of settings: Default settings, struct argument,
    %       key+value arguments.
    %       PRO: Needs to implement analogue of "add_struct_to_struct" for containers.Map.
    %   PROPOSAL: Can permit both containers.Map and struct simultaneously.
    %       PROPOSAL: Convert struct to containers.Map internally.
    
    %====================================================
    % Assign SettingsArg1: Uses first argument if struct
    %====================================================
    if numel(argList) >= 1 && isstruct(argList{1})
        SettingsArg1 = argList{1};
        argList      = argList(2:end);    % NOTE: Shortens argList.
    else
        SettingsArg1 = struct;
    end
    
    %======================================================================================
    % Assign SettingsArgListPairs: Uses remaining arguments interpreted as key-value pairs
    %======================================================================================
    SettingsArgListPairs = struct;
    while true
        if numel(argList) == 0
            break

        elseif numel(argList) == 1
            error('Uneven number of string-value arguments.')
            
        elseif numel(argList) >= 2
            if ~ischar(argList{1})
                error('Expected string argument is not string.')
            end
            
            SettingsArgListPairs.(argList{1}) = argList{2};
            
            argList = argList(3:end);
        end
    end
    
    %===============================
    % Assign return value: Settings
    %===============================
    % Take all available values from SettingsArgListPairs.
    Settings = SettingsArgListPairs;
    
    % For missing values, take from SettingsArg1.
    Settings = EJ_library.utils.add_struct_to_struct(Settings, SettingsArg1, ...
        struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'bothAreStructs', 'Recurse'));
    
    % For missing values, take from DefaultSettings.
    Settings = EJ_library.utils.add_struct_to_struct(Settings, DefaultSettings, ...
        struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'bothAreStructs', 'Recurse'));
end
