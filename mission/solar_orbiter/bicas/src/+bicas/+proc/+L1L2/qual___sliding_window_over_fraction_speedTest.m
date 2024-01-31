%
% Simple script for testing performance of function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function qual___sliding_window_over_fraction_speedTest()
close all

xAr    = [];
tSecAr = [];
%for i = logspace(log10(1e4), log10(1e4), 10)
for i = logspace(log10(1e3), log10(1e6), 10)
  nSamples           = i;
  minFlaggedFraction = 0.6;
  windowLengthSec    = 6*1;

  periodSec   = 1e99;
  samplFreqHz = 1000;

  % Log
  nSamplesPerWindow = windowLengthSec * samplFreqHz;
  fprintf('nSamplesPerWindow = %g\n', nSamplesPerWindow);

  tt2000Ar = int64( [1:nSamples]' / samplFreqHz * 1e9 );
  % square() generates a square wave with values -1 and +1 and a period of
  % 2*pi. square(0) = +1.
  bFlag1Ar = logical(square( tt2000Ar * (2*pi*periodSec/1e9)) == 1);

  if 0
    % DEBUG
    plot(double(tt2000Ar)/1e9, bFlag1Ar)
    ylim([-0.1, 1.1])
    return
  end

  xAr(end+1)    = nSamples;
  tSec          = test(tt2000Ar, bFlag1Ar, minFlaggedFraction, windowLengthSec);
  tSecAr(end+1) = tSec;

  % Log
  fprintf('------------------------\n')
  fprintf('nSamples          = %g\n', nSamples);
  fprintf('periodSec         = %g\n', periodSec);
  fprintf('tSec              = %g\n', tSec);
end

figure('WindowState','maximized')
for i = 1:2
  subplot(2, 1, i)
  if i==1
    semilogy(xAr, tSecAr, 'o-')
  elseif i==2
    loglog(  xAr, tSecAr, 'o-')
  end
  xlabel('nSamples')
  ylabel('Time [s]')
  grid on
end
end



function tSec = test(tt2000Ar, bFlag1Ar, minFlaggedFraction, windowLengthSec)
tTicToc = tic();
[~] = bicas.proc.L1L2.qual.sliding_window_over_fraction(tt2000Ar, bFlag1Ar, minFlaggedFraction, windowLengthSec);
tSec = toc(tTicToc);
end