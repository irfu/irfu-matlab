function [freqRps, z] = convert_TF_human2math(freqHz, gainEnergyDb, phaseShiftDeg)
%
% Effectively convert tabulated transfer function (TF) on conventional "human-readable" format to a "mathematically
% pure" z(omega) format.
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% freqHz        : TF frequencies.
% gainEnergyDb  : Amplitude change on the form of dB __energy__, i.e. 20 dB change corresponds to a factor
%                 of 10. Rationale: If a signal U (Volt) is amplified by X dB, then the power P~U^2
%                 (W) is amplified by 2*X dB.
% phaseShiftDeg : TF phase in degrees. (360 degrees per revolution.)
%                 NOTE: phaseShiftDeg will be uncertain up to adding n*360 degrees.
% freqRps       : TF frequencies. Unit: radians/s
% z             : Complex TF amplitudes (DFT component multiplication factors).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-12-12



freqRps = freqHz * 2*pi;   % (Convert revolutions to radians = factor 2*pi)
z       = db2mag(gainEnergyDb) .* exp(1i * deg2rad(phaseShiftDeg));

end
