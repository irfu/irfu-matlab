%
% Quick-and-dirty script for manually experimenting with
% bepic.spinfit.spin_aligned() by specifying input values and then plotting them
% and the output after processing (by editing this code).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function spin_aligned___MTEST

if 0
  % =============
  % Step function
  % =============
  tt2000Ar       = int64(20:119)';
  samplesAr      = [3*ones(1, 50)  4*ones(1, 50)]';
  spinPhaseRadAr = wrapTo2Pi(linspace(0, 20*pi, numel(tt2000Ar)))';
end

if 1
  % =========
  % Sine wave
  % =========
  SPIN_PERIOD_NS   = 4e9;
  SIGNAL_PERIOD_NS = 4e9;

  tt2000Ar         = int64(0e9 : 0.25e9 : 300e9)';
  spinPhaseRadAr   = wrapTo2Pi( double(tt2000Ar)/SPIN_PERIOD_NS * 2*pi );

  signalPhaseRadAr = double(tt2000Ar) / SIGNAL_PERIOD_NS * 2*pi;
  samplesAr        = 3 + 4*sin(signalPhaseRadAr) + 5*sin(2*signalPhaseRadAr);

  % Insert data gap on the form of NaN samples.
  %samplesAr(1000:2:1100) = NaN;
end


assert(numel(tt2000Ar) == numel(samplesAr))
assert(numel(tt2000Ar) == numel(spinPhaseRadAr))



close all
display_result( ...
  tt2000Ar            = tt2000Ar, ...
  spinPhaseRadAr      = spinPhaseRadAr, ...
  samplesAr           = samplesAr, ...
  timeWindowPeriodRad = 2*pi, ...
  timeWindowLengthRad = 2*pi, ...
  timeWindowCenterRad = 1*pi, ...
  dataGapMinNs        = int64(2e9));
end



function display_result(A)
  arguments
    A.tt2000Ar
    A.spinPhaseRadAr
    A.samplesAr
    A.timeWindowPeriodRad
    A.timeWindowLengthRad
    A.timeWindowCenterRad
    A.dataGapMinNs
  end

R = bepic.spinfit.spin_aligned(...
  tt2000Ar            = A.tt2000Ar, ...
  spinPhaseRadAr      = A.spinPhaseRadAr, ...
  samplesAr           = A.samplesAr, ...
  timeWindowPeriodRad = A.timeWindowPeriodRad, ...
  timeWindowLengthRad = A.timeWindowLengthRad, ...
  timeWindowCenterRad = A.timeWindowCenterRad, ...
  dataGapMinNs        = A.dataGapMinNs);

% R.tt2000Ar
% R.offsetAr

% ====
% PLOT
% ====
figure
t = tiledlayout(7, 1, "Padding", "compact", "TileSpacing", "compact");
ax1 = add_tile("samplesAr",         A.tt2000Ar, A.samplesAr);
ax2 = add_tile("spinPhaseRadAr",    A.tt2000Ar, A.spinPhaseRadAr);
ax3 = add_tile("offsetAr",          R.tt2000Ar, R.offsetAr);
ax4 = add_tile("coefficientCos1Ar", R.tt2000Ar, R.coefficientCos1Ar);
ax5 = add_tile("coefficientSin1Ar", R.tt2000Ar, R.coefficientSin1Ar);
ax6 = add_tile("coefficientCos2Ar", R.tt2000Ar, R.coefficientCos2Ar);
ax7 = add_tile("coefficientSin2Ar", R.tt2000Ar, R.coefficientSin2Ar);

% NOTE: Important to link X axes, since they are not identical otherwise, and
% could be decieving.
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7], 'x');
end



function ax = add_tile(titleStr, tt2000Ar, yAr)
LINE_TYPE_CA = {"-o", "LineWidth", 2};

ax = nexttile;
plot(ax, tt2000Ar, yAr, LINE_TYPE_CA{:})
grid on
title(ax, titleStr)
end