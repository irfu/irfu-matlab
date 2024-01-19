%
% Generic general-purpose function for applying a TF (linear frequency-dependent
% transfer function) to a sequence of real-valued (time domain) samples.
%
% Wrapper around bicas.tf.apply_TF_kernel().
%
%
% ARGUMENTS
% =========
% dt
%       Scalar, numeric. Seconds. Time between samples.
% y1
%       Column vector. Signal.
% tf
%       Function handle. Z = tf(omegaRps). Transfer function.
% lenKernel
%       Length of kernel that shall be generated from tf.
% edgePolicy
% varargin
%       Optional settings arguments as interpreted by
%       irf.utils.interpret_settings_args().
%
%
% RETURN VALUES
% =============
% y2
%       Column vector. Transformed signal.
% yKernelB
%       Column vector.
%       Kernel used internally, after the Hann window was applied to it (if
%       enabled).
%       NOTE: Only returned for testing purposes.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-08-21.
%
function [y2, yKernelB] = apply_TF_time(dt, y1, tf, lenKernel, edgePolicy, varargin)
    %
    % PROPOSAL: edgePolicy as setting, not named argument.
    %
    % PROPOSAL: MTEST code.
    %   PROPOSAL: Use data from L1R datasets. Compare freq. and time domain
    %           applications.
    %   PROPOSAL: Good enough plots for meetings.

    EMID = 'BICAS:Assertion:IllegalArgument';

    DEFAULT_SETTINGS.detrendingDegreeOf = -1;
    DEFAULT_SETTINGS.retrendingEnabled  = false;
    DEFAULT_SETTINGS.hannWindow         = false;

    Settings = irf.utils.interpret_settings_args(...
        DEFAULT_SETTINGS, varargin);
    irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
    clear DEFAULT_SETTINGS



    %=============
    % ~ASSERTIONS
    %=============
    % TODO-DEC: Which argument assertions should one bother to have?
    % bicas.tf.apply_TF_freq() and bicas.tf.apply_TF_kernel() check most
    % arguments.
    if ~isnumeric(lenKernel)
        error(EMID, 'lenKernel is not numeric..')
    elseif ~isscalar(lenKernel)
        error(EMID, 'lenKernel is not scalar.')
    elseif ~(lenKernel>0)
        error(EMID, 'lenKernel is not positive.')
    end
    assert(islogical(Settings.hannWindow))



    %===============
    % Obtain kernel
    %===============
    % NOTE: Set "kernel origin" (see bicas.tf.apply_TF_kernel()) to middle, and
    % rounded down for even-length kernels. The algorithm is designed so that
    % this can be set quite arbitrarily but in reality one probably wants to set
    % it around the middle index.
    iKernelOrigin           = floor(1 + (lenKernel-1)/2);
    yImpulse                = zeros(lenKernel, 1);
    yImpulse(iKernelOrigin) = 1;
    yKernel                 = bicas.tf.apply_TF_freq(dt, yImpulse, tf);
    % NOTE: Uses bicas.tf.apply_TF_freq(), BICAS' other main function for
    % applying transfer functions to signals, using FFT. Here it is only used
    % for obtaining an impulse response in the time domain, i.e. kernel.



    if Settings.hannWindow
        %=============================
        % Apply Hann window to kernel
        %=============================
        unshiftedHannWin = hann(lenKernel, 'periodic');

        % (Periodic) Hann window center index, i.e. where the Hann window=max=1,
        % for the initial Hann window produced by hann().
        % --
        % NOTE: Index is a half-integer for ODD-numbered-length PERIODIC Hann
        % windows (i.e. ODD-numbered-length kernels). Rounding is therefore
        % ~arbitrary.
        iInitialHannWinCenter = 1 + ceil((lenKernel-1)/2);   % Round up.

        % Disable?
        assert(...
            (mod(lenKernel, 2) == 1) || ...
            unshiftedHannWin(iInitialHannWinCenter) == 1)

        % Circularly shift Hann Window so that the Hann window max is at the
        % kernel origin.
        shiftedHannWin = circshift(...
            unshiftedHannWin, ...
            iKernelOrigin - iInitialHannWinCenter);

        yKernelB = yKernel .* shiftedHannWin;
    else
        % Do nothing.
        yKernelB = yKernel;
    end



    %================
    % Process signal
    %================
    Drt = bicas.tf.drt(Settings.detrendingDegreeOf, Settings.retrendingEnabled);
    y1b = Drt.detrend(y1);

    y2b = bicas.tf.apply_TF_kernel(y1b, yKernelB, iKernelOrigin, edgePolicy);

    % NOTE: Using freq. domain-TF for scaling.
    y2 = Drt.retrend(y2b, tf(0));

end
