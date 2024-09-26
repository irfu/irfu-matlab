%
% Interpret the contents of an (already read) config file as settings key-value
% pairs.
%
% NOTE: This codes defines a syntax for the configuration file. This must be
% compatible with what the bash wrapper script can handle to read what it needs
% itself (a subset).
%
%
% ARGUMENTS
% =========
% configFileRowList
%       Cell array of strings representing the content of a BICAS config file.
%       Permitted to set the same key twice or more. The last row with the
%       setting is used.
% settingsVsMap
%       containers.Map with all values on string form.
%       VS = values as strings.
%
%
% IMPLEMENTATION NOTE
% ===================
% Does not read the config file itself in order to make the code more easily
% testable. Does not modify a bicas.Settings object in order to be able to use
% one common code for interpreting numeric string values (reading settings from
% CLI arguments also requires converting string arguments to numeric values.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-01-24
%
function settingsVsMap = interpret_config_file(configFileRowsCa, L)

SETTINGS_KEY_REGEXP          = '[a-zA-Z0-9._-]+';
SETTINGS_VALUE_STRING_REGEXP = '[^"]*';

ASSIGNMENT_RE_LIST = {...
  SETTINGS_KEY_REGEXP, ' *= *', '"', ...
  SETTINGS_VALUE_STRING_REGEXP, '" *', '(#.*)?'};



% IMPLEMENTATION NOTE: Well-defined keyType and valueType good for (1)
% comparisons in automated tests, and (2) for restricting values.
settingsVsMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

for iRow = 1:numel(configFileRowsCa)
  row = configFileRowsCa{iRow};

  if ~isempty(regexp(row, '^#.*$', 'once'))
    % CASE: Row with only comments, beginning at first character of row!
    % Do nothing

  elseif ~isempty(regexp(row, '^ *$', 'emptymatch'))
    % CASE: Row contains only whitespace (or is empty).
    % Do nothing

  else
    % CASE: Row is a setting key assignment.
    [subStrList, ~, isPerfectMatch] = ...
      irf.str.regexp_str_parts(...
      row, ASSIGNMENT_RE_LIST, 'PERMIT_NON_MATCH');

    % ASSERTION
    if ~isPerfectMatch
      error(...
        'BICAS:interpret_config_file:Assertion:CannotInterpretConfigFile', ...
        'Can not interpret row %i in configuration file: "%s"', ...
        iRow, row)
    end
    key      = subStrList{1};
    valueStr = subStrList{4};

    if settingsVsMap.isKey(key)

      % NOTE: Log message is annoying when running automatic testing
      % code.
      L.logf('warning', ...
        ['Settings key "%s" is assigned a second (or more) time', ...
        ' in the config file.'], ...
        key)
    end
    settingsVsMap(key) = valueStr;
  end
end

end
