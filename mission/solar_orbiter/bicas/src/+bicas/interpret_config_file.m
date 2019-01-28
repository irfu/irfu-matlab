%
% Interpret the contents of an (already read) config file as settings key-value pairs.
%
%
% ARGUMENTS
% =========
% configFileRowList : Cell array of strings representing the content of a BICAS config file.
% settingsVsMap     : containers.Map with all values on string form. VS = values as strings.
%
%
% IMPLEMENTATION NOTE: Does not read the config file itself in order to make the code more easily testable.
% Does not modify a settings object in order to be able to use one common code for interpreting numeric string values
% (reading settings from CLI arguments also requires converting string arguments to numeric values.
%
% NOTE: This codes defines a syntax for the configuration file. This must be compatible with what the bash wrapper
% script can handle to read what it needs.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-24
%
function settingsVsMap = interpret_config_file(configFileRowList)

SETTINGS_KEY_REGEXP = '[a-zA-Z0-9._-]+';
SETTINGS_VALUE_STRING_REGEXP = '[^"]*';

% IMPLEMENTATION NOTE: Well-defined keyType and valueType good for comparisons in automated tests, and (2) for
% restricting values.
settingsVsMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

for iRow = 1:numel(configFileRowList)
    row = configFileRowList{iRow};
    
    if ~isempty(regexp(row, '^#.*$'))   % Row with only comments, beginning at first character!
        ;
    elseif ~isempty(regexp(row, '^ *$', 'emptymatch'))   % Row containing only whitespace (or empty).
        ;
    else
        try
            subStrList = bicas.utils.regexp_str_parts(row, {SETTINGS_KEY_REGEXP, ' *= *', '"', SETTINGS_VALUE_STRING_REGEXP, '" *', '(#.*)?'});
            key      = subStrList{1};
            valueStr = subStrList{4};
            
        catch exception
            error('BICAS:interpret_config_file:Assertion:CannotInterpretConfigFile', 'Can not interpret row %i in configuration file: "%s"', iRow, row)
        end
        
        if settingsVsMap.isKey(key)
            error('BICAS:interpret_config_file:Assertion:CannotInterpretConfigFile', 'The same settings key "%s" is assigned twice.', key)
        end
        settingsVsMap(key) = valueStr;
    end
end

end



