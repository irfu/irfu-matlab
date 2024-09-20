%
% Utility function for interpreting settings argument(s) on a standardized
% syntax/format passed to functions.
%
%
% INTENDED USAGE
% ==============
% The function interprets a 1D cell array, representing a list of arguments to
% another function, typically the external function that calls this function,
% and typically the external function's last arguments by using varargin.
% --
% IMPORTANT NOTE: The function does permit that settings (fields) not present in
% DefaultSettings are ADDED to Settings. In order to require/permit a certain
% set of settings (keys), use an assertion for the set of field names in the
% returned struct. Note that if there are default values for all keys, then
% DefaultSettings can be used to obtain the list of field names.
%   Ex: irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
% --
% NOTE: One could conceivably use MATLAB's inputParser instead to achieve
% similar ends, but that was not known by the autor at the time of writing this
% code. irf.utils.interpret_settings_args() also provides functionality for
% overriding settings in bulk and in sequence (submit a struct) which
% inputParser does not.
%
%
% ALGORITHM
% =========
% The algorithm uses/constructs three analogous "Settings" structs.
% DefaultSettings      = argument.
% SettingsArg1         = argsCa{1}, if argsCa{1} exists and is a struct.
%                        Otherwise empty struct.
% SettingsArgListPairs = remainder of argsCa (excluding SettingsArg1),
%                        interpreted as pairs of field name + field value.
% --
% The functions returns a struct which is a combination (union of fields) of
% (1) SettingsArgListPairs
% (2) SettingsArg1
% (3) DefaultSettings
% where field values are taken from the first top-most struct with that field,
% i.e. a higher one has precedence over a lower one, e.g. (1) has precedence
% over (2). Therefore, the structs do not need to, but are allowed to, contain
% the same fields.
%
%
% LIMITATIONS
% ===========
% Can not handle situations where set of fields in the "full settings struct"
% depends on other fields.
%   Ex: Settings.enabled = 1 implies there should be a field Settings.parameter,
%       but Settings.enabled = 0 implies there should not.
%
%
% ARGUMENTS
% =========
% DefaultSettings
%       Struct with default Settings to be processed (not Settings for this
%       function).
% argsCa
%       Cell array representing a sequence of arguments (presumably "varargin"
%       or subset thereof) from another function that uses this function.
%       It is either
%           (1) {              key1,value1, ..., keyN,valueN}, or
%           (2) {SettingsArg1, key1,value1, ..., keyN,valueN}
%       NOTE: The cell array is permitted to be empty.
%
%
% RETURN VALUES
% =============
% Settings
%       Struct. See algorithm.
% argsCa
%       Cell array of strings, representing list of arguments passed to other
%       function. Typically varargin as received from the enclosing function
%       directly.
%
%
% Initially created 2018-07-18 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function Settings = interpret_settings_args(DefaultSettings, argsCa)
% PROPOSAL: Assert SettingsArg1 and SettingsArgListPairs fields to always exist in DefaultSettings.
%   CON: Makes it impossible to have default values for some settings/fields, but not for others.
%   PROPOSAL: Option/flag for this behaviour.
%   NOTE: User can easily add an assertion after:
%       irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
%
% PROPOSAL: Reorg algorithm to produce better error messages when using bad combinations of arguments.
%
% PROPOSAL: Argument for required settings fields (DefaultSettings contains the ones which are optional in varargin).
%   CON: Might have situations where the required settings fields depend on the the value of other settings fields.
%   PROPOSAL: User should call irf.assert.struct instead. This is almost as succint.
%
% PROPOSAL: Add ability to recognize "string keywords" (one argument, instead of keyword+value) which indicate that
% a flag (false/true) shall be set.
%   NOTE: The default value is always false.
%   TODO-DEC: How represent in a struct that a field represents an optional string keyword?
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
%   PROPOSAL: Can support both containers.Map and struct simultaneously. Return same type as DefaultSettings
%             (for backward compatibility).
%       PROPOSAL: Convert struct to containers.Map internally.
%
% PROPOSAL: Policy argument (first) for whether or not to only permit settings/fields already in DefaultSettings.
%   CON: If wants to permit, then still undetermined (not asserted) which extra settings should be permitted?
%       PROPOSAL: Cell array of settings that are allowed to be added. No extra settings. ==> Empty cell array.
%   PRO: Prevents caller from forgetting this option.
%
% PROBLEM: Mixing interface (key strings) with implementation (variable/field names). Can not change fieldname
%          without changing interface.
%
% PROPOSAL: Extend syntax: Permit struct at "any" and multiple points in list of arguments.
%   Keys present in later location overwrite earlier.
%   Ex: struct1,key1,value1,struct2,key2,value2
%   CON: Has not actually had the need for this. Is just "clever".
%   NOTE: Can be implemented recursively?!!
%
% PROPOSAL: Somehow make compatible with mixing arguments to multiple
%           destinations.
%   Ex; Mix arguments for conventional settings, and plotting properties.
%   CON: Better to group arguments together in cell arrays.
%   CON: Setting name collisions.
%
% PROPOSAL: Shorter name.
%   PRO: Frequenctly used and referenced function.
%   ~parse=p
%   ~interpret=i
%   psa, isa, isargs, parse_sargs, psargs

assert(iscell(argsCa), 'Argument "argsCa" is not a cell array.')

%====================================================
% Assign SettingsArg1: Uses first argument if struct
%====================================================
if numel(argsCa) >= 1 && isstruct(argsCa{1})
  SettingsArg1 = argsCa{1};
  argsCa       = argsCa(2:end);    % NOTE: Shortens argsCa.
else
  SettingsArg1 = struct;
end

%======================================================================
% Assign SettingsArgListPairs: Uses remaining arguments interpreted as
% key-value pairs
%======================================================================
SettingsArgListPairs = struct;
while true
  if numel(argsCa) == 0
    break

  elseif numel(argsCa) == 1
    error('Uneven number of string-value arguments.')

  elseif numel(argsCa) >= 2
    if ~ischar(argsCa{1})
      error('Expected string argument is not string.')
    end

    SettingsArgListPairs.(argsCa{1}) = argsCa{2};

    argsCa = argsCa(3:end);
  end
end

%===============================
% Assign return value: Settings
%===============================
% Take all available values from SettingsArgListPairs.
%Settings = SettingsArgListPairs;

% NOTE: Could permit recursive settings structs, but then without checking
% for the existence of corresponding default key.
AS2S_SETTINGS = struct(...
  'noStructs',    'Overwrite', ...
  'aIsStruct',    'Error', ...
  'bIsStruct',    'Error', ...
  'abAreStructs', 'Error', ...
  'onlyBField',   'Copy');

% For missing values, take from SettingsArg1.
%Settings = add_struct_to_struct(Settings, SettingsArg1, AS2S_SETTINGS);
Settings = add_struct_to_struct(DefaultSettings, SettingsArg1, AS2S_SETTINGS);

% For missing values, take from DefaultSettings.
%Settings = add_struct_to_struct(Settings, DefaultSettings, AS2S_SETTINGS);
Settings = add_struct_to_struct(Settings, SettingsArgListPairs, AS2S_SETTINGS);
end



%===============================================================================
% Lightly modified hard-coded copy of
% irf.ds.add_struct_to_struct()
% -------------------------------------------------------------------------
% Modifications: Not use irf.utils.interpret_settings_args
% (including recursive call).
%
% IMPORTANT IMPLEMENTATION NOTE
% =============================
% This function (irf.utils.interpret_settings_args) DELIBERATELY
% DOES NOT USE irf.ds.add_struct_to_struct() so that it (that
% function) can use irf.utils.interpret_settings_args() INSTEAD.
% Both using each other would lead to infinite recursion.
%===============================================================================
function A = add_struct_to_struct(A, B, Settings)

%===========================================================================
% IMPLEMENTATION NOTE: CAN NOT USE
% irf.utils.interpret_settings_args() SINCE IT USES THIS
% FUNCTION!
%===========================================================================
if nargin == 2
  Settings = DEFAULT_SETTINGS;
elseif nargin == 3
  ;   % Do nothing. Settings has already been assigned.
else
  % ASSERTION
  error('Wrong number of arguments.')
end

% ASSERTIONS
assert(isstruct(A),               'A is not a structure.')
assert(isstruct(B),               'B is not a structure.')
assert(isstruct(Settings),        'Settings is not a struct.')
assert(isscalar(A),               'A is not scalar.')
assert(isscalar(B),               'B is not scalar.')
%assert(isequal(size(A), size(B)), 'A and B have different array sizes.')

%===========================
% Iterate over fields in B.
%===========================
bFieldNamesList = fieldnames(B);
for i = 1:length(bFieldNamesList)
  fieldName = bFieldNamesList{i};

  if isfield(A, fieldName)
    %===========================================================
    % CASE: Duplicate fields (field exists in both structures).
    %===========================================================
    afv = A.(fieldName);
    bfv = B.(fieldName);

    abAreStructs = false;
    if ~isstruct(afv) && ~isstruct(bfv)
      behaviour = Settings.noStructs;

    elseif isstruct(afv) && ~isstruct(bfv)
      behaviour = Settings.aIsStruct;

    elseif ~isstruct(afv) && isstruct(bfv)
      behaviour = Settings.bIsStruct;

    else
      behaviour = Settings.abAreStructs;
      abAreStructs = true;
    end

    if strcmp(behaviour, 'Error')
      error('Structures share identically named fields "%s".', fieldName)
    elseif strcmp(behaviour, 'Overwrite')
      A.(fieldName) = B.(fieldName);
    elseif strcmp(behaviour, 'Do nothing')
      % Do nothing.
    elseif strcmp(behaviour, 'Recurse') && (abAreStructs)
      % NOTE: RECURSIVE CALL. Needs the original Settings.
      A.(fieldName) = add_struct_to_struct(A.(fieldName), B.(fieldName), Settings);
    else
      error(...
        ['Can not interpret string value behaviour="%s"', ...
        ' for this combination of field values.'], ...
        behaviour)
    end
  else
    %===========================================
    % CASE: Field exists in "B", but not in "A"
    %===========================================
    switch(Settings.onlyBField)
      case 'Copy'
        % NOTE: Not overwrite field, but create new field in A.
        A.(fieldName) = B.(fieldName);
      case 'Error'
        error('Field B.%s exists but not A.%s.', fieldName, fieldName)
      otherwise
        error('Illegal setting onlyBField=%s', Settings.onlyBField)
    end

  end
end
end
