%
% Speed test code. To estimate the speed of execution.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function spinfit_fit_SAFW_pot___speedTest
S = generate_signal_sine_wave();

tic
R = bepic.spinfit.fit_SAFW_pot(...
  tt2000Ar              = S.tt2000Ar, ...
  spinPhaseRadAr        = S.spinPhaseRadAr, ...
  samplesAr             = S.samplesAr, ...
  fitWindowPeriodRad    = 2*pi, ...
  fitWindowLengthRad    = 2*pi, ...
  fitWindowCenterRefRad = 1*pi, ...
  dataGapMinNs          = int64(2e9));
toc

R;
end



% Generate example input signal.
function S = generate_signal_sine_wave()
SPIN_PERIOD_NS   = 4e9;
SIGNAL_PERIOD_NS = 4e9;
SAMPLING_RATE_HZ = 32;    % Highest sampling rate?!

tt2000Ar         = int64([0 : 1/SAMPLING_RATE_HZ : 24*60*60] * 1e9)';    % 24 h
spinPhaseRadAr   = wrapTo2Pi( double(tt2000Ar)/SPIN_PERIOD_NS * 2*pi );

signalPhaseRadAr = double(tt2000Ar) / SIGNAL_PERIOD_NS * 2*pi;

samplesAr(:, 1)  = 8 + 1*sin(signalPhaseRadAr) + 2*sin(2*signalPhaseRadAr);
samplesAr(:, 2)  = 9 + 3*cos(signalPhaseRadAr) + 4*cos(2*signalPhaseRadAr);

S = struct();
S.tt2000Ar       = tt2000Ar;
S.samplesAr      = samplesAr;
S.spinPhaseRadAr = spinPhaseRadAr;
end
