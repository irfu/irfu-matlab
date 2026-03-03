%
% Utility function for easily deriving the bias current at specified timestamps
% almost directly from SOLO_L1_RPW-BIA-CURRENT dataset zVariables. Therefore
% has some extra features to handle the quirks of the dataset format.
%
%
% NOTE: "CURRENT" in the function name refers to SOLO_L1_RPW-BIA-CURRENT CDFs.
% NOTE: Is used by BICAS, which needs to detect a known anomaly (multiple
% instances of equivalent successive bias settings) and react to it itself
% (choose between e.g. error and mitigation).
%
%
% ARGUMENTS
% =========
% t1
%       Nx1 vector. Time. Sorted. Exact type of time is unimportant as long
%       as it is a number that increases linear with time). Can be TT2000.
% zvIBIASx1
%       Nx1 vector. Floating-point. Same length as t1. New bias values.
%       Can be zVar from CURRENT dataset.
% t2
%       Nx1 vector. Time. Same time format as t1.
%
%
% RETURN VALUES
% =============
% zvIBIASx2
%       Current values at t2. Floating-point. NaN for timestamps before first
%       input timestamp.
% duplicatesAnomaly
%       Scalar logical. Whether the function detected a known anomaly in
%       SOLO_L1_RPW-BIA-CURRENT datasets.
%       RATIONALE: Caller (e.g. BICAS) can give error, warning.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
function [zvIBIASx2, duplicatesAnomaly] = CURRENT_ZV_to_current_interpolate(t1, zvIBIASx1, t2)

[t1b, zvIBIASx1b, duplicatesAnomaly] = solo.hwzv.CURRENT_ZV_to_current(t1, zvIBIASx1);



zvIBIASx2 = interpolate(t1b, zvIBIASx1b, t2);
end



% Interpolate using the nearest previous value.
%
% IMPLEMENTATION NOTE: Bias currents are set VERY RARELY. Must therefore
% interpolate by looking at previous value.
%
%

function y2 = interpolate(x1, y1, x2)
% IMPLEMENTATION NOTE: MATLAB's interp1() (with option "previous") has many
% strange requirements. Therefore not using it.
% NOTE: interp1() does NOT require t1 to be sorted.
% NOTE: interp1() requires x1 values to be unique, even if the corresponding
% y1 values are identical.
% NOTE: interp1() requires x_data to be finite (y_data can be non-finite though).
% NOTE: interp1() requires at least two data points, also when using option
% "previous".
% NOTE: MATLAB typecasts to the type of the array when indexing.
%     zvIBIASx2 = interp1(t1b, zvIBIASx1b, t2, 'previous');
% NOTE: interp1(... 'previous') returns NaN also for t2 > max(t1)!! Must
% therefore do this oneself.
%     zvIBIASx2(t2 > t1b(end)) = zvIBIASx1b(end);

% ASSERTIONS
assert(iscolumn(x1), 'x1 is not a column vector.')
assert(iscolumn(y1), 'y1 is not a column vector.')
% x2 should not need to have any particular size.
assert(issorted(x1, 'strictascend'), 'Time stamps are not increasing monotonically.')
assert(all(isfinite(x1)))



% NOTE: x1 must be sorted for algorithm to work.

% Pre-allocate. Set default value for timestamps before first timestamp.
y2 = NaN(size(x2));   % NOTE: Always double.

if numel(x1) >= 1
  for i = 1:numel(x1)-1
    b     = (x1(i) <= x2) & (x2 < x1(i+1));
    y2(b) = y1(i);
  end

  b     = x2 >= x1(end);
  y2(b) = y1(end);
end
end
