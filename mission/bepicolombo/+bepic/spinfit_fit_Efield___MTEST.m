%
% Script for manually experimenting with fitting E field by specifying input
% values and then plotting them and the output after processing (by editing this
% code).
%
% This code is intended for testing multiple such E field functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function spinfit_fit_Efield___MTEST
% PROPOSAL: See as usable for multiple forms of fitting?

%R = generate_signal_step_function();
%R = generate_signal_sine_wave();
R = generate_signal_from_file();

assert(numel(R.tt2000Ar) == numel(R.samplesAr))
assert(numel(R.tt2000Ar) == numel(R.spinPhaseRadAr))

close all
if 1
  fit_display_result___fit_SAFW_MMS( ...
    tt2000Ar              = R.tt2000Ar, ...
    spinPhaseRadAr        = R.spinPhaseRadAr, ...
    samplesAr             = R.samplesAr, ...
    fitWindowPeriodRad    = 2*pi*1, ...
    fitWindowLengthRad    = 2*pi*1, ...
    fitWindowCenterRefRad = 2*pi*0.5, ...
    nMinFitSamples        = 6, ...
    dataGapMinNs          = int64(2e9));
end
if 0
  fit_display_result___fit_SAFW_E( ...
    tt2000Ar              = R.tt2000Ar, ...
    spinPhaseRadAr        = R.spinPhaseRadAr, ...
    samplesAr             = R.samplesAr, ...
    fitWindowPeriodRad    = 2*pi*1, ...
    fitWindowLengthRad    = 2*pi*1, ...
    fitWindowCenterRefRad = 2*pi*0.5, ...
    dataGapMinNs          = int64(2e9));
end
end



% ###############
% EXAMPLE SIGNALS
% ###############



% Generate example input signal.
function R = generate_signal_step_function()
R.tt2000Ar       = int64(20:119)';
R.samplesAr      = [3*ones(1, 50)  4*ones(1, 50)]';
R.spinPhaseRadAr = wrapTo2Pi(linspace(0, 20*pi, numel(tt2000Ar)))';
end



% Generate example input signal.
function R = generate_signal_sine_wave()
SPIN_PERIOD_NS   = 4e9;
SIGNAL_PERIOD_NS = 4e9;

tt2000Ar         = int64([0 : 0.25 : 100] * 1e9)';

% Create data gap on the form of MISSING INDICES (not NaN).
b = (20e9 < tt2000Ar) & (tt2000Ar < 30e9);
tt2000Ar = tt2000Ar(~b);

spinPhaseRadAr   = wrapTo2Pi( double(tt2000Ar) / SPIN_PERIOD_NS * 2*pi );

signalPhaseRadAr = double(tt2000Ar) / SIGNAL_PERIOD_NS * 2*pi;
samplesAr        = 3 + 4*sin(signalPhaseRadAr) + 2*cos(2*signalPhaseRadAr);

% Create data gap on the form of NaN SAMPLES.
b = (70e9 < tt2000Ar) & (tt2000Ar < 80e9);
samplesAr(b) = NaN;

R = struct();
R.tt2000Ar       = tt2000Ar;
R.samplesAr      = samplesAr;
R.spinPhaseRadAr = spinPhaseRadAr;
end



% Generate example input signal.
function R = generate_signal_from_file()
FILE = '/media/erjo/bepicolombo/reformatter_mirror/pwi/cdf/EFD/L1_prime/2025/bc_mmo_pwi-efd_l1p_l-e_20251016_r01-v00-00.cdf';
Do = dataobj(FILE);

tt2000Ar       = bepic.spinfit.utils.E_ZV_epoch_to_linear(...
  Do.data.epoch.data, ...
  Do.data.t_offset_4hz.data);
spinPhaseRadAr = bepic.spinfit.utils.E_ZV_spinphase_to_linear(...
  Do.data.spinphase_4hz.data);
samplesAr      = bepic.spinfit.utils.E_ZV_Ev_to_linear(...
  Do.data.Ev_4hz.data);

R.tt2000Ar       = tt2000Ar;
R.samplesAr      = samplesAr;
R.spinPhaseRadAr = spinPhaseRadAr;
end



% ##########
% FIT & PLOT
% ##########



% Function for calculating the spin fit and plotting some of the spin fit input
% and output.
function fit_display_result___fit_SAFW_MMS(A)
  arguments
    A.tt2000Ar
    A.spinPhaseRadAr
    A.samplesAr
    A.fitWindowPeriodRad
    A.fitWindowLengthRad
    A.fitWindowCenterRefRad
    A.nMinFitSamples
    A.dataGapMinNs
  end

tic
R = bepic.spinfit.fit_SAFW_MMS(...
  tt2000Ar              = A.tt2000Ar, ...
  spinPhaseRadAr        = A.spinPhaseRadAr, ...
  samplesAr             = A.samplesAr, ...
  fitWindowPeriodRad    = A.fitWindowPeriodRad, ...
  fitWindowLengthRad    = A.fitWindowLengthRad, ...
  fitWindowCenterRefRad = A.fitWindowCenterRefRad, ...
  nMinFitSamples        = A.nMinFitSamples, ...
  dataGapMinNs          = A.dataGapMinNs, ...
  nFitCoefficients      = 5);
toc

% ====
% PLOT
% ====
figure
t = tiledlayout(7, 1, "Padding", "compact", "TileSpacing", "compact");

axAr = matlab.graphics.axis.Axes.empty(0, 1);
%axAr(end+1) = add_tile("INPUT: tt2000Ar",               A.tt2000Ar, A.tt2000Ar);
axAr(end+1) = add_tile("INPUT: samplesAr",              A.tt2000Ar, A.samplesAr);
axAr(end+1) = add_tile("INPUT: spinPhaseRadAr",         A.tt2000Ar, A.spinPhaseRadAr);
%axAr(end+1) = add_tile("OUTPUT: fitWindowCenterTt2000", R.fitWindowCenterTt2000, R.fitWindowCenterTt2000);
axAr(end+1) = add_tile("OUTPUT: offsetAr",          R.fitWindowCenterTt2000, R.offset);
axAr(end+1) = add_tile("OUTPUT: coefficientCos1Ar", R.fitWindowCenterTt2000, R.coefficientCos1);
axAr(end+1) = add_tile("OUTPUT: coefficientSin1Ar", R.fitWindowCenterTt2000, R.coefficientSin1);
axAr(end+1) = add_tile("OUTPUT: coefficientCos2Ar", R.fitWindowCenterTt2000, R.coefficientCos2);
axAr(end+1) = add_tile("OUTPUT: coefficientSin2Ar", R.fitWindowCenterTt2000, R.coefficientSin2);

% NOTE: Important to link X axes, since they are not identical otherwise, and
% could be deceiving.
linkaxes(axAr(:)', 'x');
end



% Function for calculating the spin fit and plotting some of the spin fit input
% and output.
function fit_display_result___fit_SAFW_E(A)
  arguments
    A.tt2000Ar
    A.spinPhaseRadAr
    A.samplesAr
    A.fitWindowPeriodRad
    A.fitWindowLengthRad
    A.fitWindowCenterRefRad
    A.dataGapMinNs
  end

tic
R = bepic.spinfit.fit_SAFW_E(...
  tt2000Ar              = A.tt2000Ar, ...
  spinPhaseRadAr        = A.spinPhaseRadAr, ...
  samplesAr             = A.samplesAr, ...
  fitWindowPeriodRad    = A.fitWindowPeriodRad, ...
  fitWindowLengthRad    = A.fitWindowLengthRad, ...
  fitWindowCenterRefRad = A.fitWindowCenterRefRad, ...
  dataGapMinNs          = A.dataGapMinNs);
toc

% ====
% PLOT
% ====
figure
t = tiledlayout(7, 1, "Padding", "compact", "TileSpacing", "compact");

axAr = matlab.graphics.axis.Axes.empty(0, 1);
%axAr(end+1) = add_tile("INPUT: tt2000Ar",               A.tt2000Ar, A.tt2000Ar);
axAr(end+1) = add_tile("INPUT: samplesAr",              A.tt2000Ar, A.samplesAr);
axAr(end+1) = add_tile("INPUT: spinPhaseRadAr",         A.tt2000Ar, A.spinPhaseRadAr);
%axAr(end+1) = add_tile("OUTPUT: fitWindowCenterTt2000", R.fitWindowCenterTt2000, R.fitWindowCenterTt2000);
axAr(end+1) = add_tile("OUTPUT: A(1)=offsetAr",          R.fitWindowCenterTt2000, R.A(:, 1));
axAr(end+1) = add_tile("OUTPUT: A(2)=coefficientCos1Ar", R.fitWindowCenterTt2000, R.A(:, 2));
axAr(end+1) = add_tile("OUTPUT: A(3)=coefficientSin1Ar", R.fitWindowCenterTt2000, R.A(:, 3));
axAr(end+1) = add_tile("OUTPUT: A(4)=coefficientCos2Ar", R.fitWindowCenterTt2000, R.A(:, 4));
axAr(end+1) = add_tile("OUTPUT: A(5)=coefficientSin2Ar", R.fitWindowCenterTt2000, R.A(:, 5));

% NOTE: Important to link X axes, since they are not identical otherwise, and
% could be deceiving.
linkaxes(axAr(:)', 'x');
end



% ##################
% REUSABLE PLOT CODE
% ##################



function ax = add_tile(titleStr, tt2000Ar, yAr)
LINE_TYPE_CA = {"-o", "LineWidth", 2};

ax = nexttile;
plot(ax, tt2000Ar, yAr, LINE_TYPE_CA{:})
grid on
title(ax, titleStr)
end
