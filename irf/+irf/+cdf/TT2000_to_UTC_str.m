%
% Convert one CDF TT2000 value to one UTC string with nanoseconds resolution.
%
% NOTE: spdfparsett2000() seems to work well as an inverse without wrapper.
%
%
% ARGUMENT
% ========
% tt2000
%       Scalar numeric value. Does not have to be int64.
%
%
% RETURN VALUE
% ============
% utcStr
%       String (char array). Example: 2020-04-01T01:23:45.678901234
% nSecondDecimals
%       Number of decimals with which the number of seconds should be printed.
%       Printed value is rounded (not truncated).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-03.
%
function utcStr = TT2000_to_UTC_str(tt2000, nSecondDecimals)
% TODO-DEC: How handle various needs for formats? Rounding, truncation?
% PROPOSAL: Assertions on argument being int64 as they are in CDF?
% NOTE: Should be analogous to any inverted conversion function.

assert(isscalar(tt2000), 'Illegal argument tt2000. Must be scalar.')

utcStrCa = irf.cdf.TT2000_to_UTC_str_many(tt2000, nSecondDecimals);
utcStr = utcStrCa{1};
end
