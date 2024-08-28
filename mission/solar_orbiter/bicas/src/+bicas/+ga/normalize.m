%
% Normalize an arbitrary value by replacing specified
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function x = normalize(x, beforeValuesCa, afterValue)
% PROPOSAL: Convert to generic function.
% NOTE: Compare irf.utils.translate().
%
% PROBLEM: Function is used partly  for normalizing CDF input. Can not easily
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
