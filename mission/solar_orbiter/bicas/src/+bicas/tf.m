%
% Code for applying TF in general. In practice, one main function with helper
% function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef tf
  % PROPOSAL: Automatic test code for bicas.tf.apply_TF_with_DRT()



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Apply transfer function to signal according to chosen algorithm.
    % Potentially modify TF and data before and after application of TF.
    %
    %
    % TERMINOLOGY
    % ===========
    % De-trending : REMOVING fit on data BEFORE applying the TF. It does NOT
    %               automatically imply RE-trending.
    % Re-trending : ADDING BACK a scaled version of the previously removed fit
    %               when de-trending, AFTER applying the TF.
    %
    %
    % De-trending / Re-trending
    % =========================
    % NOTE: Detrending makes it impossible to modify the amplitude & phase for
    % the frequency components in the trend itself (the fit), e.g. to delay the
    % signal in the trend. If the input signal is interpreted as N-periodic
    % (e.g. when applying TF using FFT), then de-trending affects the jump
    % between the beginning and end of the signal (reduces it in the case of
    % linear or higher order de-trending), which should reduce erroneous
    % high-frequency content (counterexample: sine wave, one cycle). The
    % implementation scales the "trend" (polynomial fit) by tfZ(omega==0).
    % --
    % NOTE: Retrending is bad for non-lowpass filters (TFs) since the
    % retrending requires scaling the fit by tfZ(omega=0) which is only
    % meaningful for lowpass filters.
    % --
    % ** Code has the ability to enable/disable de-trending:
    %       -- To handle both DC and AC signals.
    %       -- Make testing easier.
    %
    %
    % ARGUMENTS
    % =========
    % dtSec
    %       Scalar. Time difference between samples in seconds.
    % y1
    %       Column vector of samples. May be modified by this function before
    %       actually applying the TF.
    % tf
    %       Function handle. Transfer function. Same as for
    %       bicas.tf.freq.apply_TF() and bicas.tf.time.apply_TF().
    % S
    %       Optional settings arguments. Available settings:
    %         * 'detrendingDegreeOf'
    %               Integer.
    %               >=0 : Degree of the polynomical fit used for de-trending.
    %                <0 : No de-trending.
    %               Default = -1.
    %         * 'retrendingEnabled'
    %         * 'snfSubseqMinSamples'
    %               Positive integer. Required minimum number of finite-valued
    %               samples in a subsequence to processed (not returned as NaN).
    %         +more. See implementation.
    %
    %
    % RETURN VALUES
    % =============
    % y2
    %       y1 after applying modifications and TF. Always column vector.
    % Debug
    %       Struct. For debugging and automatic tests.
    %       .y1ModifCa
    %           Potentially modified y1 on which the TF is applied.
    %       .y2ModifCa
    %           Data after applying TF and before it is potentially modified
    %           (again).
    %       .tfModif
    %           The actual TF used. "Modif" in the name is for historical
    %           reasons.
    %       .i1Array, .i2Array
    %           1D column arrays with indices into y1 for the first and last
    %           index for the respective time intervals separated by non-finite
    %           y1 values.
    %
    %
    % Author: Erik P G Johansson, Uppsala, Sweden
    % First created 2020-11-04.
    %
    function [y2, Debug] = apply_TF(dtSec, y1, tf, S)
      % PROPOSAL: Move bicas.tf to bicas.proc.L1L2.tf.
      %   PRO: Code is only used for processing L1/L1R-->L2.
      %   CON: Implies that code is less generic/reusable.
      %
      % PROPOSAL: Function for splitting up samples separated by non-finite
      %           samples (SNF). -- IMPLEMENTED
      % PROPOSAL: This function should not support SNF. Should be made easy by
      %           separate function instead.
      %   PROBLEM: Needs way of easily putting subsequences together again after
      %            calibration.
      %     PROPOSAL: Function which (given algorithm) splits array into CA of
      %               arrays representing *ALL* samples, including ones which
      %               can not be calibrated.
      %               One can then easily execute below steps:
      %               (1) Use function to split array-->CA of arrays.
      %               (2) Iterate over arrays and calibrate each separately and
      %                   store the result in new CA of arrays.
      %                   NOTE: It is assumed that subsequences which can not be
      %                   calibrated (due to containing invalid values) are
      %                   calibrated and that the calibration code fails
      %                   nicely (e.g. producing only NaN).
      %               (3) Combine calibrated CA of arrays into one array.
      %               -- IMPLEMENTED
      %       PRO: Can be implemented using irf.utils.split_by_change().
      %       CON: Calibrates NaN sequences --> NaN. Multiplies the number of
      %            sequences to calibrate by ~2.
      %
      % PROPOSAL: Check that data is finite. Only call bicas.tf.freq.apply_TF()
      %           if all data is non-finite. -- IMPLEMENTED
      %   PRO: bicas.tf.freq.apply_TF() can assume (needs to be updated) that
      %        always Z<>NaN and thereby detect if TF can not be evaluated via NaN.
      %       PRO: Can construct TFs in steps/parts where each part does not have
      %            to be evaluated for all omega (return NaN if can not be
      %            evaluated).
      %           CON: Not necessarily best solution. TFs could give error when
      %                not being able to return value.
      %
      % PROPOSAL: Always snfEnabled=true. Remove setting.
      %
      % PROPOSAL: Better abbbreviation for SNF.
      %   ~split, values, samples, time series
      %   FV
      %     Reflects the sample values in the actual datasets, i.e. what the user
      %     sees. Relevant if using term for settings.
      %     FVs do not represent +-inf.
      %   Non-finite
      %     Reflects the sample values in the input to this function, what the
      %     function actually sees.
      %   Invalid values
      %     Keeps the criterion for splitting open.
      %   --
      %   SNF  = Split by Non-Finite -- IMPLEMENTED
      %   SNFS = Split by Non-Finite Samples
      %   SFV = Split by FV
      %   NFS = Non-Finite Splitting
      %   SBIV = Split by Invalid Values
      %
      % PROPOSAL: Remove Debug.
      %   NOTE: Only used in bicas.tf___UTEST.

      arguments
        dtSec, y1, tf
        S.detrendingDegreeOf      = -1;
        S.retrendingEnabled       = false;
        S.method                  = 'FFT';
        S.kernelLengthMax         = Inf;
        S.kernelEdgePolicy        = 'MIRROR';
        S.kernelHannWindowEnabled = false;
        S.snfEnabled              = false;
        S.snfSubseqMinSamples     = 1;
      end

      assert(...
        isscalar( S.snfEnabled) && ...
        islogical(S.snfEnabled))
      assert(...
        isscalar(S.snfSubseqMinSamples) && ...
        S.snfSubseqMinSamples >= 1)

      % ASSERTION: Arguments
      assert(isscalar(dtSec) & (dtSec > 0))
      assert(iscolumn(y1), 'Argument y1 is not a column vector.')



      if S.snfEnabled
        %===================================================================
        % Split up time interval into sub-intervals separated by non-finite
        % samples (fill values)
        %===================================================================
        y1Ca = bicas.tf.split_samples_by_nonfinite(y1);
      else
        y1Ca = {y1};
      end

      % Pre-allocate and initialize values that will not later be overwritten.
      % y2  = NaN(size(y1));

      nSs = numel(y1Ca);

      Debug = struct();
      Debug.y1ModifCa = cell(nSs, 1);   % Pre-allocate
      Debug.y2ModifCa = cell(nSs, 1);   % Pre-allocate
      Debug.tfModif   = tf;

      y2Ca = cell(nSs, 1);
      for iSs = 1:nSs
        y1Ss = y1Ca{iSs};
        if numel(y1Ss) < S.snfSubseqMinSamples
          y1Ss = NaN(size(y1Ss));
        end
        [y2ss, D] = bicas.tf.apply_TF_with_DRT(dtSec, y1Ss, tf, S);

        Debug.y1ModifCa{iSs} = D.y1Modif;
        Debug.y2ModifCa{iSs} = D.y2Modif;

        y2Ca{iSs} = y2ss;
      end
      y2 = cell2mat(y2Ca);
      y2 = y2(:);   % Normalize 0x0 --> 0x1.

    end    % function apply_TF()



    function yCa = split_samples_by_nonfinite(y)
      % NOTE: Does not implement any constraint on the minimum length of
      %       subsequences since that can easily be implemented separately.
      %       Also does not want to implement special behaviour for such
      %       subsequences, e.g. set to NaN.

      % PROPOSAL: Better function name
      %   not "samples"?
      %     CON: Is not a true, generic function. Does not return indices (like
      %          irf.utils.split_by_change()), but a CA of arrays.
      %   split_by_nonfinite()
      %     CON: Sounds too generic.
      %
      % PROPOSAL: Permit arbitrary "criterion", bArray.

      assert(iscolumn(y) & isnumeric(y))

      [i1Array, i2Array] = irf.utils.split_by_change(isfinite(y));

      n   = length(i1Array);
      yCa = cell(n, 1);
      for i = 1:n
        yCa{i, 1} = y(i1Array(i):i2Array(i));
      end

    end



    % Apply TF using detrending.
    function [y2, Debug] = apply_TF_with_DRT(dt, y1, tf, Settings)
      %-------------------------------------------------
      % Ensure that data is either all valid or all NaN
      %-------------------------------------------------
      % NOTE: Also ensures identical behaviour for FFT and KERNEL methods wrt.
      % to this (though that is probably not a problem anyway).
      if any(~isfinite(y1))
        y1 = NaN(size(y1));
      end



      %#####################
      % Optionally DE-trend
      %#####################
      Drt = bicas.tf.Deretrending(...
        Settings.detrendingDegreeOf, ...
        Settings.retrendingEnabled);
      y1Detrended = Drt.detrend(y1);



      %#########################
      % APPLY TRANSFER FUNCTION
      %#########################
      switch(Settings.method)

        case 'FFT'
          y2Detrended = bicas.tf.freq.apply_TF(dt, y1Detrended, tf);

        case 'KERNEL'
          % Capping kernel length to improve performance.
          lenKernel = min(length(y1), Settings.kernelLengthMax);

          % NOTE: The called function applies the Hann window instead of current
          % function since it only applies to kernel method (as opposed to
          % de-trending & re-trending).
          y2Detrended = bicas.tf.time.apply_TF(...
            dt, y1Detrended, tf, lenKernel, Settings.kernelEdgePolicy, ...
            Settings.kernelHannWindowEnabled);

        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            'Illegal setting "method" value.')
      end



      %#####################
      % Optionally RE-trend
      %#####################
      y2 = Drt.retrend(y2Detrended, tf(0));



      Debug = struct();
      Debug.y1Modif = y1Detrended;
      Debug.y2Modif = y2Detrended;
    end



    % Create modified TF where Z=0 above certain frequency.
    %
    % ARGUMENTS
    % =========
    % tfHighFreqLimitFraction
    %       Fraction of Nyquist frequency (1/dt). TF is set to zero above this
    %       frequency. Can be Inf.
    %
    function tf2 = make_hard_low_pass_TF(tf, tfHighFreqLimitFraction, dtSec)
      % PROPOSAL: Move function to bicas.tf.utest_utils and redefine that file
      %           as bicas.tf.utils.
      %   NOTE: This function is *NOT* called from inside bicas.tf (e.g.
      %         bicas.tf.apply_TF()).

      % NOTE: Permits tfHighFreqLimitFraction to be +Inf.
      assert(...
        isnumeric(  tfHighFreqLimitFraction) ...
        && isscalar(tfHighFreqLimitFraction) ...
        && ~isnan(  tfHighFreqLimitFraction) ...
        && (        tfHighFreqLimitFraction >= 0))

      % Nyquist frequency [rad/s] =
      % = 2*pi [rad/sample] * (1/2 * 1/dt [samples/s])
      % = pi/dt
      nyquistFreqRps     = pi/dtSec;
      tfHighFreqLimitRps = tfHighFreqLimitFraction * nyquistFreqRps;
      tf2                = ...
        @(omegaRps) (tf(omegaRps) .* (omegaRps < tfHighFreqLimitRps));
    end



  end    % methods(Static)



end
