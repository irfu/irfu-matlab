%
% Speed test code. To estimate the speed of execution.
%
%
% CONCLUSIONS
% ===========
% 2026-11-20, EJ:
%   sarcasm, 24 h, 32 Hz ==> Execution takes 17 s.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function spinfit_fit_pot_SAFW___speedTest
format compact

% ===============
% Run speed tests
% ===============
lengthSecAr = logspace(log10(0.1*60*60), log10(24*60*60), 10);
wallTimeSecAr = [];
for lengthSec = lengthSecAr
  lengthSec
  wallTimeSecAr(end+1) = get_time_consumption(lengthSec)
end
% Log the combined results
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
title("bepic.spinfit.fit_pot_SAFW()", "interpreter", "none")
end



function wallTimeSec = get_time_consumption(lengthSec)
S = generate_signal_sine_wave(lengthSec);

tic
R = bepic.spinfit.fit_pot_SAFW(...
  tt2000Ar                = S.tt2000Ar, ...
  spinPhaseRadAr          = S.spinPhaseRadAr, ...
  samplesAr               = S.samplesAr, ...
  dataGapMinNs            = int64(2e9), ...
  nMinSamplesPerFitWindow = 1);
wallTimeSec = toc;
end



% Generate example input signal.
function S = generate_signal_sine_wave(lengthSec)
SPIN_PERIOD_NS   = 4e9;
SIGNAL_PERIOD_NS = 4e9;
SAMPLING_RATE_HZ = 32;    % Highest sampling rate for potential data?!

tt2000Ar         = int64([0 : 1/SAMPLING_RATE_HZ : lengthSec] * 1e9)';    % 24 h
spinPhaseRadAr   = wrapTo2Pi( double(tt2000Ar)/SPIN_PERIOD_NS * 2*pi );

signalPhaseRadAr = double(tt2000Ar) / SIGNAL_PERIOD_NS * 2*pi;

samplesAr(:, 1)  = 8 + 1*sin(signalPhaseRadAr) + 2*sin(2*signalPhaseRadAr);
samplesAr(:, 2)  = 9 + 3*cos(signalPhaseRadAr) + 4*cos(2*signalPhaseRadAr);

S = struct();
S.tt2000Ar       = tt2000Ar;
S.samplesAr      = samplesAr;
S.spinPhaseRadAr = spinPhaseRadAr;
end
