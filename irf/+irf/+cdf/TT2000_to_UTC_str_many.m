%
% Convert CDF tt2000 value to UTC string with nanoseconds.
%
% NOTE: spdfparsett2000 seems to work well as an inverse without wrapper.
%
%
% ARGUMENT
% ========
% tt2000 : Arbitrary size array of numeric values. Does not have to be int64.
%
%
% RETURN VALUE
% ============
% utcStrCa : Cell array of strings. Same size as tt2000. Example: '2020-04-01T01:23:45.678901234'
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-03.
%
function utcStrCa = TT2000_to_UTC_str_many(tt2000Array)
% TODO-DEC: How handle various needs for formats? Rounding, truncation?
% PROPOSAL: Assertions on argument being int64 as they are in CDF?
% PROPOSAL: Handle array.
%   NOTE: Return value must be cell array.
%   NOTE: Special case for empty array unless backward-incompatible.
%   NOTE: Useful to be able to return string, not cell. See actually made calls.
%   PROPOSAL: Separate function TT2000_to_UTC_str_many.
%       CON: Might duplicate future functionality.
%       PROPOSAL: Implement this function using *many function.
%   PROPOSAL: Flag for returning cell or string.
%
% NOTE: Should be analogous to any inverted conversion function.

% NOTE: irf.cdf.TT2000_to_datevec can handle Nx1 arrays, where N>=1.
dateVec = irf.cdf.TT2000_to_datevec(tt2000Array(:));

utcStrCa = cell(size(tt2000Array));
for i = 1:numel(utcStrCa)
  utcStrCa{i} = sprintf('%04i-%02i-%02iT%02i:%02i:%012.9f', dateVec(i, 1:6));
end
end
