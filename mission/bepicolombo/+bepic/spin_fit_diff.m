%
% Function which implements the core processing for MEFISTO spin fits, when
% processing L1p diff CDFs --> L2pre diff CDFs for constant-length fit time
% windows.
%
% The purpose of this file is to isolate many of the technical details
% surrounding the details of the spin fit, e.g. the time windows used for spin
% fitting, and which exact timestamps should represent spin fit time windows.
%
% This function is essentially a wrapper around mms.mms_spinfit_m. In the long
% term, this function should probably be replaced by a copy dedicated to
% BepiColombo processing (in case it needs to be modified).
%
%
% TECHNICALITIES WHICH ARE BEYOND THE SCOPE OF THIS FUNCTION
% ==========================================================
% * Rotating the coordinate system.
% * Converting first-order sin/cos fit coefficients to a vector.
%
%
% ARGUMENTS
% =========
% tt2000Ar
%       TT2000 timestamps for samples.
% spinPhaseRadiansAr
%       Spin phase values
% samplesAr
%       Sample values.
% --
% All arguments are same-length column arrays.
%
%
% RETURN VALUES
% =============
% Struct with self-explanatory fieldnames. One spin fit per row.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function R = spin_fit_diff(A)
% PROPOSAL: Better name
%   ~diff, E field
%   MEFISTO, Mio
%   Type of time windows: constant-length, spin phase-aligned
%   Does this function really "know" that the samples represent a "diff"?
% PROPOSAL: Put in package dedicated to official processing, bepic.*.
% PROPOSAL: Put in class with other methods for processing.
%   Ex: Spin fit for constant-length and spin-phase aligned time windows.
%   Ex: Processing diff and potential data.
%   PRO: Can share constants and ~private functions.
% PROPOSAL: Separate wrapper for calling mms_spinfit_m().
%   PRO: Handles special case of zero-length arrays.
%   PRO: Enforces type-checking.
%   PRO: Names return value coefficients
%     CON: The number of return value coefficients depends on number of terms
%          (input).
%       CON: Unlikely to change?
%
% PROPOSAL: Support both constant-length and spin phase-aligned time windows.
%
% PROPOSAL: Expose constants are arguments?
%   Ex: Length of time windows.
%   CON: Need to write more tests.
%     PRO: Might need to test for behaviour which is not used.
%   PRO: Easier to change constants (update implementation)
%
% TODO-NI: Spin phase in radians?
% TODO-DEC: Range of spin phase? Infinite? 0 to 2*pi.
%   PROPOSAL: Assert spin phase interval: 0 to 2*pi; -pi to +pi
%     PRO: Checks radians vs degrees
% 
%
% PROPOSAL: Automatic test code.
%
% [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 5, int64([1:6]'), sin([1:6]'), [1:6]', 4e9, 4e9, 0)
% n=10; spinRad=linspace(0, 0.1*pi, n)'; [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 3+2, int64([1:n]'), 3+cos(spinRad)+2*sin(spinRad), spinRad, 4e9, 4e9, 0)
% n=10; spinRad=linspace(0, 0.1*pi, n)'; [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 3+2, int64([1:n]'), 3+cos(spinRad)+2*sin(spinRad), spinRad, 4e9, 4e9, 0)

arguments
  A.tt2000Ar
  A.spinPhaseRadiansAr
  A.samplesAr
end

N_FIT_TERMS                        = 3+2;
N_MAX_FIT_ITERATIONS               = 5;              % TODO: Determine proper value.
N_MIN_REQUIRED_FIT_SAMPLES         = N_FIT_TERMS+3;  % TODO: Determine proper value.
NS_BETWEEN_TIME_WINDOWS_BEGINS     = 4e9;            % Nanoseconds
TIME_WINDOW_LENGTH_NS              = NS_BETWEEN_TIME_WINDOWS_BEGINS;
TIME_WINDOW_BEGIN_REFERENCE_TT2000 = int64(0);



tt2000Ar           = A.tt2000Ar;
spinPhaseRadiansAr = A.spinPhaseRadiansAr;
samplesAr          = A.samplesAr;

assert(iscolumn(tt2000Ar)           & isa(tt2000Ar,           'int64'))
assert(iscolumn(spinPhaseRadiansAr) & isa(spinPhaseRadiansAr, 'double') & all(isfinite(spinPhaseRadiansAr)))
assert(iscolumn(samplesAr)          & isa(samplesAr,          'double'))
nIn = numel(tt2000Ar);
assert(nIn == numel(spinPhaseRadiansAr))
assert(nIn == numel(samplesAr))



% ===============================================
% COPIED FROM THE mms.mms_spinfit_m DOCUMENTATION
% ===============================================
% function [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(maxIt, minPts, nTerms, timeData, data, phase, fitEvery, fitInterv, t0)
%  Compute spinfit coefficients to spinning data. Data is fitted to
%  function y = A + Bcos(phase) + Csin(phase) + (Dcos(2*phase) +
%  Esin(2*phase) + Fcos(3*phase) + Gsin(3*phase)). According to the number
%  of terms specified.
%
% Input: (all required)
%   maxIt    - maximum of iterations for each fit
%   minPts   - minimum number of points required for each fit
%   nTerms   - number of terms to fit, must be odd (3, 5, 7)
%   timeData - time of measurement (int64 TT2000)
%   data     - data (samples) to be fitted
%   phase    - phase of instrument at corresponding time of measurement
%              (radians)
%   fitEvery - one spinfit every X ns (default every 5*10^9 ns)
%   fitInter - spinfit is fitted to data during this interval
%              (default 20*10^9 ns)
%   t0       - the first time inside timeData which is evenly divisable
%              with fitEvery, accounting for leap seconds and such.
%              (With default, each fit should line up with times 00:00:05,
%              00:00:10, 00:00:15 etc), (int64 TT2000).
% Output: (all required)
%   timeFit  - middle of each spinfit ( 00:00:05, 00:00:10 etc). (int64 TT2000)
%   sfit     - matrix with each fit coefficents
%   sdev     - standard deviation of each fit
%   iter     - number of iterations used for each fit
%   nBad     - number of bad points, outliers for each fit
%
% Bad fits will have value NaN.
%
% This is an interface function used by Matlab to display help and/or
% hints, the real processing occurs in mms_spinfit_mx (mex file).
% ==============================================================================
if nIn == 0
  % ========================
  % CASE: inpuy empty arrays
  % ========================
  % IMPLEMENTATION NOTE: Using empty arrays when calling mms_spinfit_m() can
  % crash MATLAB!!! ('25.2.0.3042426 (R2025b) Update 1', Linux). Can therefore
  % not call it for this case.  
  % Ex: [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 5, int64.empty(0, 1), double.empty(0, 1), double.empty(0, 1), 4e9, 4e9, 0)

  timeFit = int64.empty( 0, 1);
  sfit    = double.empty(0, 5);
  sdev    = double.empty(0, 1);
  iter    = double.empty(0, 1);
  nBad    = double.empty(0, 1);
else
  % ============================
  % CASE: Non-empty input arrays
  % ============================
  [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(...
    N_MAX_FIT_ITERATIONS, N_MIN_REQUIRED_FIT_SAMPLES, N_FIT_TERMS, ...
    tt2000Ar, samplesAr, spinPhaseRadiansAr, ...
    NS_BETWEEN_TIME_WINDOWS_BEGINS, TIME_WINDOW_LENGTH_NS, TIME_WINDOW_BEGIN_REFERENCE_TT2000);
end

% =========================================
% Assertions on mms_spinfit_m return values
% =========================================
% IMPLEMENTATION NOTE: Useful for (1) testing the understanding/behaviour of
% mms_spinfit_m(), and (2) ensuring that emulated mms_spinfit_m() return values
% are consistent with the non-emulated ones.
nOut = numel(timeFit);
assert(iscolumn(timeFit) & (numel(timeFit) == nOut))
assert(isequal(size(sfit), [nOut, N_FIT_TERMS]))
assert(iscolumn(sdev)    & (numel(sdev)    == nOut))
assert(iscolumn(iter)    & (numel(iter)    == nOut))
assert(iscolumn(nBad)    & (numel(nBad)    == nOut))



% =========================
% Construct return value(s)
% =========================
R = struct();
R.offsetAr          = sfit(:, 1);
R.coefficientCos1Ar = sfit(:, 2);
R.coefficientSin1Ar = sfit(:, 3);
R.coefficientCos2Ar = sfit(:, 4);
R.coefficientSin2Ar = sfit(:, 5);
end