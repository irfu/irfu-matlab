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



    % Apply transfer function to signal according to chosen algorithm. Potentially
    % modify TF and data before and after application of TF.
    %
    %
    % TERMINOLOGY
    % ===========
    % De-trending : REMOVING fit on data BEFORE applying the TF. It does NOT
    %               automatically imply RE-trending.
    % Re-trending : ADDING BACK a scaled version of the previously removed fit when
    %               de-trending, AFTER applying the TF.
    %
    %
    % De-trending / Re-trending
    % =========================
    % NOTE: Detrending makes it impossible to modify the amplitude & phase for the
    % frequency components in the trend itself (the fit), e.g. to delay the signal
    % in the trend. If the input signal is interpreted as N-periodic (e.g. when
    % applying TF using FFT), then de-trending affects the jump between the
    % beginning and end of the signal (reduces it in the case of linear or higher
    % order de-trending), which should reduce erroneous high-frequency content
    % (counterexample: sine wave, one cycle). The implementation scales the "trend"
    % (polynomial fit) by tfZ(omega==0).
    % --
    % NOTE: Retrending is bad for non-lowpass filters (TFs) since the retrending
    % requires scaling the fit by tfZ(omega=0) which is only meaningful for lowpass
    % filters.
    % --
    % ** Code has the ability to enable/disable de-trending:
    %       -- To handle both DC and AC signals.
    %       -- Make testing easier.
    % ** Code has the ability to make TF zero above frequency cutoff. This cut-off
    %    is naturally sampling frequency-dependent and is therefore NOT a natural
    %    part of the TF itself.
    %
    %
    % ARGUMENTS
    % =========
    % dt
    %       Scalar. Time difference between samples.
    % y1
    %       Column vector of samples. May be modified by this function before
    %       actually applying the TF.
    % tf
    %       Function handle. Transfer function. Same as for bicas.tf.freq.apply_TF().
    % S
    %       Optional settings arguments. Available settings:
    %         * 'detrendingDegreeOf'
    %               Integer.
    %               >=0 : Degree of the polynomical fit used for de-trending.
    %                <0 : No de-trending.
    %               Default = -1.
    %         * 'retrendingEnabled'
    %         * 'tfHighFreqLimitFraction'
    %               Fraction of Nyquist frequency (1/dt). TF is regarded as zero
    %               above this frequency. Can be Inf.
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
    %           The actual (potentially modified) TF used.
    %       .i1Array, .i2Array
    %           1D column arrays with indices into y1 for the first and last index
    %           for the respective time intervals separated by non-finite y1 values.
    %
    %
    % Author: Erik P G Johansson, Uppsala, Sweden
    % First created 2020-11-04.
    %
    function [y2, Debug] = apply_TF(dt, y1, tf, S)
      % PROPOSAL: Move bicas.tf to bicas.proc.L1L2.tf.
      %   PRO: Code is only used for processing L1/L1R-->L2.
      %   CON: Implies that code is less generic/reusable.
      %
      % PROPOSAL: Function for splitting up samples separated by non-finite
      %           samples (SNF).
      %
      % PROPOSAL: Check that data is finite. Only call bicas.tf.freq.apply_TF()
      %           if all data is non-finite.
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
      %   --
      %   SNF  = Split by Non-Finite -- IMPLEMENTED
      %   SNFS = Split by Non-Finite Samples
      %   SFV = Split by FV
      %   NFS = Non-Finite Splitting

      arguments
        dt, y1, tf
        S.detrendingDegreeOf      = -1;
        S.retrendingEnabled       = false;
        S.tfHighFreqLimitFraction = Inf;
        S.method                  = 'FFT';
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
      assert(isscalar(dt) & (dt > 0))
      assert(iscolumn(y1), 'Argument y1 is not a column vector.')



      %=========================================================================
      % Create modified version of TF which is set to zero for high frequencies
      %=========================================================================
      % NOTE: Permits Settings.tfHighFreqLimitFraction to be +Inf.
      assert(...
        isnumeric(  S.tfHighFreqLimitFraction) ...
        && isscalar(S.tfHighFreqLimitFraction) ...
        && ~isnan(  S.tfHighFreqLimitFraction) ...
        && (        S.tfHighFreqLimitFraction >= 0))
      % Nyquist frequency [rad/s] =
      % = 2*pi [rad/sample] * (1/2 * 1/dt [samples/s])
      % = pi/dt
      nyquistFreqRps     = pi/dt;
      tfHighFreqLimitRps = S.tfHighFreqLimitFraction * nyquistFreqRps;
      tfModif            = ...
        @(omegaRps) (tf(omegaRps) .* (omegaRps < tfHighFreqLimitRps));



      if S.snfEnabled
        %===================================================================
        % Split up time interval into sub-intervals separated by non-finite
        % samples (fill values)
        %===================================================================
        % SS = SubSequence
        [i1Array, i2Array] = irf.utils.split_by_false(isfinite(y1));
        nSs                = numel(i1Array);
      else
        i1Array = 1;
        i2Array = numel(y1);
        nSs     = 1;
      end

      % Pre-allocate and initialize values that will not later be overwritten.
      y2 = NaN(size(y1));

      Debug = struct();
      Debug.y1ModifCa = cell(nSs, 1);   % Pre-allocate
      Debug.y2ModifCa = cell(nSs, 1);   % Pre-allocate
      Debug.i1Array   = i1Array;
      Debug.i2Array   = i2Array;
      Debug.tfModif   = tfModif;

      for iSs = 1:nSs
        i1 = i1Array(iSs);
        i2 = i2Array(iSs);

        y1ss = y1(i1:i2);

        if numel(y1ss) >= S.snfSubseqMinSamples
          [y2ss, D] = bicas.tf.apply_TF_with_DRT(dt, y1ss, tfModif, S);

          Debug.y1ModifCa{iSs} = D.y1Modif;
          Debug.y2ModifCa{iSs} = D.y2Modif;
        else
          y2ss = NaN(size(y1ss));

          Debug.y1ModifCa{iSs} = [];
          Debug.y2ModifCa{iSs} = [];
        end

        y2(i1:i2) = y2ss;
      end

    end    % function apply_TF()



    % Apply TF using detrending.
    function [y2, Debug] = apply_TF_with_DRT(dt, y1, tf, Settings)

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
          % PROBLEM: Kernel length == Signal length
          %          ==> Bad performance for very long time series.
          % NOTE: Length also affects amount of allocated memory (kernel, padding).
          lenKernel = length(y1);
          %lenKernel = min(lenKernel, lenKernelMax);

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



  end    % methods(Static)



end
