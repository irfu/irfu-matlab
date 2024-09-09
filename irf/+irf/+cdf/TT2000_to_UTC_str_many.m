%
% Convert CDF TT2000 value to UTC string.
%
% NOTE: spdfparsett2000() seems to work well as an inverse without wrapper.
%
%
% ARGUMENT
% ========
% tt2000Array
%       Arbitrary size array of numeric values. Does not have to be int64.
%
%
% RETURN VALUE
% ============
% utcStrCa
%       Cell array of strings (char arrays). Same array size as argument
%       "tt2000".
%       Example: '2020-04-01T01:23:45.678901234Z'
% nSecondDecimals
%       Number of decimals with which the number of seconds should be printed.
%       Printed value is rounded (not truncated).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-03.
%
function utcStrCa = TT2000_to_UTC_str_many(tt2000Array, nSecondDecimals)
% PROPOSAL: Assertions on argument being int64 as they are in CDF?
%   NOTE: bicas.utils.TT2000_to_UTC_str() adds such an assertion.
%   CON: Not needed. The value only needs to be a number, conceptually. Meaning
%        is clear without int64.
%
% NOTE: Should be analogous to any future inverted conversion function.

assert(isscalar(nSecondDecimals) && isnumeric(nSecondDecimals) ...
  && nSecondDecimals >= 0)

% NOTE: irf.cdf.TT2000_to_datevec() can only handle Nx1 arrays, where N>=1.
dateVec = irf.cdf.TT2000_to_datevec(tt2000Array(:));

% Beginning of sprintf() pattern which is always used.
PATTERN_YMDHM_PREFIX = '%04i-%02i-%02iT%02i:%02i:';

if nSecondDecimals == 0
  pattern = [PATTERN_YMDHM_PREFIX, '%02.0fZ'];
else
  nChars = 3 + nSecondDecimals;
  pattern = [PATTERN_YMDHM_PREFIX, ...
    '%0', num2str(nChars), '.', num2str(nSecondDecimals), 'fZ'];
end

utcStrCa = cell(size(tt2000Array));
for i = 1:numel(utcStrCa)
  utcStrCa{i} = sprintf(pattern, dateVec(i, 1:6));
end

end
