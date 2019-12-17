% str = sprint_SETTINGS()
%
% Create human-readable multi-line string to represent SETTINGS. Meant for logging and printing to stdout.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
function str = sprint_SETTINGS(SETTINGS)

% PROPOSAL: Somehow print where values come from (default, CLI, config file).
% PROPOSAL: Make hierarchy visually clearer?!!! Should then have help from data structure itself.

% IMPLEMENTATION NOTE: Only prints "Settings" as a header (not "constants") to indicate/hint that it is only the content
% of the "SETTINGS" variables, and not of constants.m.
str =       sprintf('\nSETTINGS\n');
str = [str, sprintf(  '========\n')];

keyList = sort(SETTINGS.get_keys());   % Values seem sorted from the method, but sort again just to be sure.
lengthMaxKey = max(cellfun(@length, keyList));



for iKey = 1:length(keyList)
    key   = keyList{iKey};
    value = SETTINGS.get_fv(key);
    
    if ischar(value)
        strValue = ['"', value, '"'];
    elseif isnumeric(value)
        strValue = sprintf('%d ', value);    % Extra whitespace important for printing arrays. Works for all dimensionalities (are made into "string row vector").
    end
    
    
    %isDefaultValue = SETTINGS.is_default_value(key);
    valueSource = SETTINGS.get_value_source(key);
    valueStatusStr = EJ_library.utils.translate({...
        {'default'},            '  --';
        {'configuration file'}, '(conf)'; 
        {'CLI arguments'},      '(CLI)'}, ...
        valueSource, 'BICAS:sprintf_settings:Assertion', 'Illegal setting value source');
    
    str = [str, sprintf(['%-6s  %-', int2str(lengthMaxKey),'s = %s\n'], valueStatusStr, key, strValue)];
end

str = [str, newline];
str = [str, sprintf('Explanations for leftmost column above:\n')];
str = [str, sprintf('---------------------------------------\n')];
str = [str, sprintf('  --   = Default value\n')];
str = [str, sprintf('(conf) = Value comes from configuration file\n')];
str = [str, sprintf('(CLI)  = Value comes from CLI argument\n')];

str = [str, newline];

end
