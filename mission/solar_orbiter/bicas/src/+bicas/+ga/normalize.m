%
% Given a value, if it is identical to any one in a specified list (column cell
% array) of values, replace it with a specified value.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function x = normalize(x, beforeValuesCa, afterValue)
% PROPOSAL: Rename to not clash with other normalization.
%   Ex: Adding "none" to empty list of strings.
% PROPOSAL: Move to new class bicas.ga (keep path).
% PROPOSAL: Convert to generic function.
% NOTE: Compare irf.utils.translate().
%
% PROBLEM: Function is used partly for normalizing CDF input. Can not easily
% distinguish between non-compliant CDFs (and give warning/error) and
% normalization.

assert(...
  iscell(beforeValuesCa) & iscolumn(beforeValuesCa), ...
  'beforeValuesCa is not a column cell array.')

for i = 1:numel(beforeValuesCa)
  beforeValue = beforeValuesCa{i};

  if isequaln(x, beforeValue)
    x = afterValue;
    return
  end
end

end
