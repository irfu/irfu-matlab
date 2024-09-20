%
% Convert CDF TT2000 values to MATLAB date vectors (UTC).
%
% Effectively a wrapper around spdfbreakdowntt2000() (TT2000-->Nx9).
%
%
% ARGUMENTS
% =========
% tt2000
%       Nx1 numeric array. N>=0. Does not have to be int64.
%
%
% RETURN VALUE
% ============
% dateVecUtc
%       Nx6 numeric array. (i,:) = MATLAB date vector. See MATLAB function
%       "datevec".
%       NOTE: (:,6) are typically decimal numbers and therefore double.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-25.
%
function dateVecUtc = TT2000_to_datevec(tt2000)
% PROPOSAL: UTC in function name.
% PROPOSAL: Separate wrapper function to handle special cases for
%           spdfbreakdowntt2000() rather than when converting to datevec
%           specifically.

assert(iscolumn(tt2000), 'Argument tt2000 is not a column array.')

% CASE: tt2000 is a Nx1 array.

if ~isempty(tt2000)
  % NOTE: spdfbreakdowntt2000() requires at least one row.
  % NOTE: spdfbreakdowntt2000() accepts both double and int64.

  % OUT = spdfbreakdowntt2000(tt2000) returns the UTC date/time from CDF
  % TT2000 time. OUT is an array with each row having nine (9) numerical
  % values for year, month, day, hour, minute, second, millisecond,
  % microsecond and nanosecond.
  v = spdfbreakdowntt2000(tt2000);
else
  v = zeros(0,9);
end
% CASE: v is a Nx9 array.

% IMPLEMENTATION NOTE: Rounding to nine digits is good for automatic
% testing. Eliminates some of arbitrariness in double-precision calculation.
v6         = round(v(:,6) + 1e-3*v(:,7) + 1e-6*v(:,8) + 1e-9*v(:,9), 9);
dateVecUtc = [v(:,1:5), v6];
end
