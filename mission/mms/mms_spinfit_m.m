function [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(maxIt, minPts, nTerms, timeData, data, phase, fitEvery, fitInterv, t0)
%  Compute spinfit coefficients to spinning data. Data is fitted to
%  function y = A + Bcos(phase) + Csin(phase) + (Dcos(2*phase) +
%  Esin(2*phase) + Fcos(3*phase) + Gsin(3*phase)). According to the number
%  of terms specified.
% Input: (all required)
%   maxIt    - maximum of iterations for each fit
%   minPts   - minimum number of points required for each fit
%   nTerms   - number of terms to fit, must be odd (3, 5, 7)
%   timeData - time of measurement (int64 TT2000)
%   data     - data to be fitted
%   phase    - phase of instrument at corresponding time of measurement
%   fitEvery - one spinfit every X ns (default every 5*10^9 ns)
%   fitInter - spinfit is fitted to data during this interval (default 20*10^9 ns)
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

narginchk(9,9);
nargoutchk(5,5);

% Ensure input is in the proper format, (double and columns).
timeData = double(timeData(:)');
data = double(data(:)');
phase = double(phase(:)');
fitEvery = double(fitEvery);
fitInterv = double(fitInterv);
t0 = double(t0);

% Call the mex function.
[timeFit, sfit, sdev, iter, nBad] = mms_spinfit_mx(maxIt, minPts, nTerms, ...
  timeData, data, phase, fitEvery, fitInterv, t0);

% Replace FillValue -159e7 with proper NaN
sfit(sfit==-159e7) = NaN;
sdev(sdev==-159e7) = NaN;
iter(iter==-159e7) = NaN;
nBad(nBad==-159e7) = NaN;

% Flip it to row
timeFit = int64(timeFit(:)); % int64 TT2000 times
sfit = sfit';
sdev = sdev(:);
iter = iter(:);
nBad = nBad(:);

end