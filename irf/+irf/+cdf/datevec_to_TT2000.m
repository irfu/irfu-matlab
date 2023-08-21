%
% Convert MATLAB date vector (UTC) to CDF TT2000.
%
% Effectively a wrapper around spdfcomputett2000 (Nx9--TT2000).
%
%
% ARGUMENTS
% =========
% dateVecUtc : Nx6 numeric array. N>=0. (i,:) = MATLAB date vector. See MATLAB
%              function "datevec".
%
%
% RETURN VALUE
% ============
% tt2000 : Nx1 array of int64.
%          Empirically: Can be special value (int64(-inf)+3) if
%          spdfcomputett2000 fails to convert, e.g. for out-of-range date.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-xx.
%
function tt2000 = datevec_to_TT2000(dateVecUtc)
% PROPOSAL: Separate wrapper function to handle spdfcomputett2000 special case.

% PROPOSAL: Make handle failed conversions when spdfcomputett2000 returns a too low value.
%   ~CON: spdfcomputett2000 already returns special value.

irf.assert.sizes(dateVecUtc, [NaN,6])

% IMPLEMENTATION NOTE: Using integer precision is good for automatic testing
% (empirically it is needed for predictable, reversible results).
nsod = int64(dateVecUtc(:, 6) * 1e9);    % NSOD = Nanoseconds Of Day
v6789 = irf.utils.mixed_radix.integer_to_MRD(nsod, [1000, 1000, 1000, 61]');
v6789 = fliplr(v6789);

v = [dateVecUtc(:, 1:5), v6789];

if ~isempty(v)
  %     OUT = spdfcomputett2000(datetime) returns the CDF TT2000 time. OUT is
  %     a vector of integer values of mxINT64_CLASS (int64).
  %
  %       datetime             An array with each row having nine (9) numerical
  %                            values for year, month, day, hour, minute, second,
  %                            millisecond, microsecond and nanosecond.
  %
  % IMPLEMENTATION NOTE: spdfcomputett2000 gives the WRONG return value if
  %                      using e.g. int64, int32 argument!!
  % IMPLEMENTATION NOTE: spdfcomputett2000 requires at least one row.
  tt2000 = spdfcomputett2000(double(v));
else
  tt2000 = int64(int64(ones(0,1)));
end
end
