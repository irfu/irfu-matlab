%
% Wrapper around bicas.utils.apply_TF_freq that (potentially) modifies data and
% TF.
%
%
% NOTES
% =====
% NOTE: Detrending makes it impossible to modify the amplitude & phase for the
% frequency components in the trend (the fit), e.g. to delay the signal. If the
% input signal is interpreted as N-periodic, then de-trending affects the jump
% between the beginning and end of the signal (reduces it in the case of linear
% de-trending), which affects the low-frequency content(?) but probably in a
% good way. The implementation scales the "trend" (polynomial fit) by
% tfZ(omega==0).
% --
% NOTE: Retrending is bad for non-lowpass filters since the retrending requires
% scaling the fit by tfZ(omega=0) which is only meaningful for lowpass filters.
% -- Has the ability to enable/disable de-trending to make testing easier.
% -- Has the ability to make TF zero above cutoff. This cut-off is naturally
%    sampling frequency-dependent and therefore not a natural part of the TF
%    itself.
%
%
% TERMINOLOGY
% ===========
% De-trending : REMOVING fit on data before applying the TF. It does NOT
%               automatically imply RE-trending.
% Re-trending : ADDING BACK a scaled version of the previously removed fit when
%               de-trending.
%
%
% ARGUMENTS
% =========
% dt, y1, tf
%       Same as for bicas.utils.apply_TF_freq(). May be modified by this
%       function before actually being submitted.
% varargin
%       Optional settings arguments as interpreted by
%       EJ_library.utils.interpret_settings_args.
%       Possible settings:
%         * detrendingDegreeOf
%               >=0 : Degree of the polynomical fit used for de-trending.
%               <0  : No de-trending.
%               Default = -1.
%         * retrendingEnabled
%         * tfHighFreqLimitFraction
%               Fraction of Nyquist frequency (1/dt). TF is regarded as zero
%               above this frequency. Can be Inf.
%
%
% RETURN VALUES
% =============
% y2
%       y1 after algorithm.
% y1B
%       Potentially modified y1 on which the TF is applied.
%       For debugging and automatic tests.
% y2B
%       Data after applying TF and before it is potentially modified.
%       For debugging and automatic tests.
% tfB
%       The actual (ptoentially modified) TF used.
%       For debugging and automatic tests.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-11-04.
%
function [y2, y1B, y2B, tfB] = apply_TF_freq_modif(dt, y1, tf, varargin)
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Better name.
    % NOTE: Could switch out the internal apply_TF_freq for other function, e.g. apply_TF_time.
    %
    % PROPOSAL: Return modified TF actually used.
    % PROPOSAL: Return modified y1 actually used.
    % PROPOSAL: Return struct.
    %   PRO: Avoid confusing return arguments.
    % PROPOSAL: Separate function for modifying TF.
    %   NOTE: tfHighFreqLimitFraction depends on sampling frequency and can not
    %         be done in advance.
    % PROPOSAL: One setting enableDetrendingRetrending = [enableDetrending, enableDetrending].
    %
    % PROPOSAL: Check that data is finite. Only call bicas.utils.apply_TF_freq
    %           if all data is non-finite.
    %   PRO: bicas.utils.apply_TF_freq can assume (needs to be updated) that always Z<>NaN and thereby detect if
    %   TF can not be evaluated via NaN.
    %       PRO: Can construct TFs in steps/parts where each part does not have
    %       to be evaluated for all omega (return NaN if can not be evaluated).
    %           CON: Not necessarily best solution. TFs could give error when
    %           not being able to return value.
    
    DEFAULT_SETTINGS.detrendingDegreeOf      = -1;
    DEFAULT_SETTINGS.tfHighFreqLimitFraction = Inf;
    DEFAULT_SETTINGS.retrendingEnabled       = 0;
    
    Settings = EJ_library.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
    EJ_library.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
    clear DEFAULT_SETTINGS
    
    
    
    assert(isscalar(Settings.detrendingDegreeOf))
    detrendingEnabled = (Settings.detrendingDegreeOf >= 0);
    
    assert(isscalar(Settings.retrendingEnabled))
    % Assert no RE-trending if no DE-trending.
    assert(detrendingEnabled || ~Settings.retrendingEnabled, ...
        'apply_TF_freq_modif:Assertion:IllegalArgument', ...
        'Illegal combination of settings "detrendingDegreeOf" and "retrendingEnabled".')
    
    
    
    %=========================================================================
    % Create modified version of TF which is set to zero for high frequencies
    %=========================================================================
    % NOTE: Permit Settings.tfHighFreqLimitFraction to be +Inf.
    assert(...
        isnumeric(  Settings.tfHighFreqLimitFraction) ...
        && isscalar(Settings.tfHighFreqLimitFraction) ...
        && ~isnan(  Settings.tfHighFreqLimitFraction) ...
        && (        Settings.tfHighFreqLimitFraction >= 0))
    % Nyquist frequency [rad/s] =
    % = 2*pi [rad/sample] * (1/2 * 1/dt [samples/s])
    % = pi/dt
    nyquistFreqRps     = pi/dt;
    tfHighFreqLimitRps = Settings.tfHighFreqLimitFraction * nyquistFreqRps;
    tfB = @(omegaRps) (tf(omegaRps) .* (omegaRps < tfHighFreqLimitRps));


    
    %#####################
    % Optionally DE-trend
    %#####################
    if detrendingEnabled
        nSamples = length(y1);
        trendFitsCoeffs1 = polyfit((1:nSamples)', y1, Settings.detrendingDegreeOf);
        yTrend1          = polyval(trendFitsCoeffs1, (1:nSamples)');
        y1B              = y1 - yTrend1;
    else
        y1B = y1;
    end
    
    
    
    %#########################
    % APPLY TRANSFER FUNCTION
    %#########################
    [y2B] = bicas.utils.apply_TF_freq(dt, y1B, tfB);
    %[y2B, tfOmegaLookups, tfZLookups] = bicas.utils.apply_TF_freq(dt, y1B, tfB);
    
%     % ASSERTIONS
%     if all(isfinite(y1B)) && ~all(isfinite(tfZLookups))
%         % IMPLEMENTATION NOTE: bicas.utils.apply_TF_freq (deliberately) does not
%         % throw error if TF returns non-finite values. This assertion should
%         % catch that case instead.
%         error('BICAS:apply_TF_freq_modif:Assertion', ...
%             'Failed to evaluate transfer function for all necessary frequencies.')
%     end
    
    
    
    %#####################
    % Optionally RE-trend
    %#####################
    if Settings.retrendingEnabled
        % Use Z(omega=0) to scale trend, including higher order polynomial
        % components.
        
        trendFitsCoeffs2 = trendFitsCoeffs1 * tfB(0);
        
        yTrend2 = polyval(trendFitsCoeffs2, (1:nSamples)');
        y2      = y2B + yTrend2;
    else
        y2 = y2B;
    end
    
end
