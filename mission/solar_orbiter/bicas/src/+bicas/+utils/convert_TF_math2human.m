function [freqHz, gainEnergyDb, phaseShiftDeg] = convert_TF_math2human(freqRps, z)
    %
    % Effectively convert tabulated transfer function (TF) on "mathematically pure" z(omega) format, to the conventional
    % corresponding "human-readable" quantitities.
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
    % IMPLEMENTATION NOTE
    % ===================
    % Only plotting and code for human analysis should need to prevent phaseShiftDeg from wrapping, or add n*360 degrees
    % and then in a custom-made way. Therefore this code does not, and should not, do that.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2017-12-12
    
    freqHz        = freqRps / (2*pi);
    gainEnergyDb  = mag2db(abs(z));
    phaseShiftDeg = rad2deg(angle(z));     % NOTE: "angle" returns angle in interval [-pi, pi].
    
end
