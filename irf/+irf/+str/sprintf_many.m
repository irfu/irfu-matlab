%
% Generate one separate string using sprintf for every component in numeric
% array.
%
% NOTE: Can only handle sprintf patterns with ONE variable.
%
%
% RATIONALE
% =========
% sprintf() does not appear to be able to do this. sprintf() applied to an array
% will produce one combined, ~repeating (except for value) formatted string.
% With this function, one can choose to join the separate strings with any
% arbitrary seperator, e.g. comma (without adding after the last formatted
% string).
%
%
% ARGUMENTS
% =========
% sprintfFormat
%       sprintf format string for one variable value (numeric, string), i.e.
%       containing ONE %i, or ONE %s etc.
% array
%       Numeric array or cell array, of scalar values to be sent to sprintf.
%       RATIONALE: Cell array is useful for submitting string values.
%
%
% RETURN VALUE
% ============
% stringCa : Cell array of strings. Same size as array.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-11
%
function stringCa = sprintf_many(sprintfFormat, array)
% PROPOSAL: Generalize to taking arbitrary arrays (per sprintf call), one component per value required by sprintf pattern.
%   NOTE: Can not use arrayfun/cellfun without repackaging arrays into one cell array, where every cell contains
%   cell array of values for one sprintf call.
%   TODO-DEC: Interface?
%       PROPOSAL: Interface: One N-dim cell array per sprintf pattern variable.
%           PRO: Natural extension.
%           PRO: Likely convenient for caller.
%           CON: Can not easily use with cellfun, arrayfun.
%       PROPOSAL: Interface: One N-dim cell array. Every component is in turn cell
%           array with one component per sprintf pattern variable.
%   PROPOSAL: Repackage arrays into one cell array of cell arrays.
%   PROPOSAL: Code itselfs iterate over all components.
%       CON: Difficult for arbitrary number of dimensions, if want return value with same output indices as input indices.
%           PROPOSAL: reshape+sub2ind+ind2sub.
%   PROPOSAL: Use one dimension/index of "array" for submitting multiple values to for each sprintf call.

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
