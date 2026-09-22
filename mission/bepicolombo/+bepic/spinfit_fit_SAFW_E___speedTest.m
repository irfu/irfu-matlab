%
% Speed test code for bepic.spinfit.fit_SAFW_E(). To estimate the speed of
% execution and how it scales.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function spinfit_fit_SAFW_E___speedTest

if 0
  get_time_consumption(10*60*60)
  return
end

lengthSecAr = logspace(log10(0.1*60*60), log10(24*60*60), 10);
wallTimeSecAr = [];
for lengthSec = lengthSecAr
  lengthSec
  wallTimeSecAr(end+1) = get_time_consumption(lengthSec)
end
lengthSecAr
wallTimeSecAr
wallTimeSecAr ./ lengthSecAr

% ===========
% Plot result
% ===========
close all
plot(lengthSecAr, wallTimeSecAr, "-o")
grid on
xlabel("Length of data[s]")
ylabel("Wall time [s]")
title("bepic.spinfit.fit_SAFW_E()", "interpreter", "none")
end



function wallTimeSec = get_time_consumption(lengthSec)
S = generate_signal_sine_wave(lengthSec);

%profile on

tic
R = bepic.spinfit.fit_SAFW_E(...
  tt2000Ar              = S.tt2000Ar, ...
  spinPhaseRadAr        = S.spinPhaseRadAr, ...
  samplesAr             = S.samplesAr, ...
  fitWindowPeriodRad    = 2*pi*1, ...
  fitWindowLengthRad    = 2*pi*1, ...
  fitWindowCenterRefRad = 2*pi*0.5, ...
  dataGapMinNs          = int64(2e9));
wallTimeSec = toc;

%profile off
%profile viewer
end



% Generate example input signal.
function S = generate_signal_sine_wave(lengthSec)
SPIN_PERIOD_NS   = 4e9;
SIGNAL_PERIOD_NS = 4e9;
SAMPLING_RATE_HZ = 128;    % Highest sampling rate.

%tt2000Ar         = int64([0 : 1/SAMPLING_RATE_HZ : 24*60*60] * 1e9)';    % 24 h
tt2000Ar         = int64([0 : 1/SAMPLING_RATE_HZ : lengthSec] * 1e9)';    % TEMP
spinPhaseRadAr   = wrapTo2Pi( double(tt2000Ar) / SPIN_PERIOD_NS * 2*pi );

signalPhaseRadAr = double(tt2000Ar) / SIGNAL_PERIOD_NS * 2*pi;
samplesAr        = 3 + 4*sin(signalPhaseRadAr) + 2*cos(2*signalPhaseRadAr);

S = struct();
S.tt2000Ar       = tt2000Ar;
S.samplesAr      = samplesAr;
S.spinPhaseRadAr = spinPhaseRadAr;
end
