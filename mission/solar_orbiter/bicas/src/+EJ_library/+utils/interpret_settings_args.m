%
% Function for interpreting settings argument(s) passed to functions. This function interprets a cell array of strings
% representing a list of arguments to another function, typically the external calling function, and typically the last
% argument(s) by using varargin.
% 
%
% ALGORITHM
% =========
% Works with up to three analogous "Settings" structs.
% SettingsArg1         = argList{1}, if argList{1} exists and is a struct. Otherwise empty struct.
% SettingsArgListPairs = remainder of argList (excluding SettingsArg1), interpreted as pairs of field name + field value.
% DefaultSettings      = argument.
% --
% Returns struct which is a combination of 
% (1) SettingsArgListPairs
% (2) SettingsArg1
% (3) DefaultSettings
% where fields are taken from the first top-most struct with that field, i.e. a higher one has precedence over a lower
% one, e.g. (1) has precedence over (2).
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
% argList         : Cell array representing a sequence of arguments (varargin presumably) from another
%                   function that uses this function.
%                   It is either
%                       (1) a cell array of pairs string+value, or
%                       (2) a cell array of exactly one struct, or
%                       (3) an empty cell array
%                   NOTE: The list is permitted to be empty.
%                   NOTE: {{}} is interpreted as empty list just as varargin uses it to represent no arguments.
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
    %       EJ_library.utils.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS))
    %
    % PROPOSAL: Reorg algorithm to produce better error messages when using bad combinations of arguments.
    %
    % PROPOSAL: Argument for required settings fields (DefaultSettings contains the ones which are optional in varargin).
    %   CON: Might have situations where the required settings fields depend on the the value of other settings fields.
    %   PROPOSAL: User should call EJ_library.utils.assert.struct instead. This is almost as succint.
    
    import EJ_library.*
    
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
    Settings = utils.add_struct_to_struct(Settings, SettingsArg1, ...
        struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'bothAreStructs', 'Recurse'));
    
    % For missing values, take from DefaultSettings.
    Settings = utils.add_struct_to_struct(Settings, DefaultSettings, ...
        struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'bothAreStructs', 'Recurse'));
end
