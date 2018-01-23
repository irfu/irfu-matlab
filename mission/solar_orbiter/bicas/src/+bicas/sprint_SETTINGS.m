% str = sprint_SETTINGS()
%
% Create human-readable multi-line string to represent SETTINGS. Meant for logging and printing to stdout.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
function str = sprint_SETTINGS()

% QUESTION: Should include header/title?
% PROPOSAL: Better handling of different data types (MATLAB classes). Let settings.m do conversions to strings?!
% PROPOSAL: Somehow print where values come from (default, CLI, config file).
% PROPOSAL: Automatically determine longest key length to use for column width.
% PROPOSAL: Make hierarchy visually clearer?!!! Should then have help from settings.m.

global SETTINGS

str = sprintf('\nSettings & constants:\n');

keyList = sort(SETTINGS.get_keys());   % Values seem sorted from the method, but sort again just to be sure.
for iKey = 1:length(keyList)
    key   = keyList{iKey};
    value = SETTINGS.get(key);
    
    if ischar(value)
        strValue = ['"', value, '"'];
    elseif isnumeric(value)
        strValue = sprintf('%d ', value);    % Extra whitespace important for printing arrays. Works for all dimensionalities (are made into "string row vector").
    end
    
    str = [str, sprintf('    %-60s = %s\n', key, strValue)];
end

end
