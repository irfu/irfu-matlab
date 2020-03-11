%
% Generate one separate string using sprintf for every component in numeric array.
%
% RATIONALE
% =========
% sprintf does not appear to be able to do this. sprintf applied to an array will produce one combined, ~repeating
% (except for value) formatted string. With this function, one can choose to join the separate strings with any
% arbitrary seperator, e.g. comma (without adding after the last formatted string).
%
%
% ARGUMENTS
% =========
% sprintfFormat : Format string for one variable value (numeric, string)
% array         : Numeric array or cell array, of scalar values to be sent to sprintf.
%                 RATIONALE: Cell array is useful for submitting string values.
%
%
% RETURN VALUE
% ============
% stringCa : Cell array of strings. Same size as array.
%
%
% Author: Erik P G Johansson
% First created 2020-03-11
%
function stringCa = sprintf_many(sprintfFormat, array)
    
    if iscell(array)
        stringCa = cellfun(...
            @(x) (sprintf(sprintfFormat, x)), ...
            array, 'UniformOutput', false);
        
    elseif isnumeric(array)
        
        stringCa = arrayfun(...
            @(x) (sprintf(sprintfFormat, x)), ...
            array, 'UniformOutput', false);
    end
end
