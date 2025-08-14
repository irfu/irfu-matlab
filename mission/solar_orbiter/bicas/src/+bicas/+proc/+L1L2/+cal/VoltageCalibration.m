%
% An instance of this class contains
%   (1) relevant settings (loaded from BSO) on how to calibrate voltages, and
%   (2) VCDS that can return actual calibration data (from RCTs).
% An instance may or may not contain calibration data for __ALL__ types of
% data/RCTs depending on how it was initialized.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR,
% TDS-CWF or TDS-RSWF) is identical (in all relevant parts) for both the RODP
% and ROC-SGSE pipeline.
%
%
% SHORTCOMINGS(?)
% ===============
% Does not implement parasitic capacitance due to lack of calibration values
% (at least). Should not need to implement support for this effect according to
% Thomas Chust(?).
%
%
% IMPLEMENTATION NOTES
% ====================
% Class is implemented with extra care, and therefore
% - has extra many assertions
% - is extra careful to include units in identifiers
% - is extra careful to use well defined terms/shortenings/naming conventions
% since
% (1) undiscovered calibration bugs could be considered extra bad,
% (2) it is expected to be hard to detect certain bugs,
% (3) it could be hard to use automatic testing here,
% (4) to detect changing RCT formats, in particular in RCTs from non-BIAS
%     teams.
% --
% NOTE: All calibration functions of measured data are assumed to accept data
% from all BLTSs (1-5), i.e. including TDS, in order to reduce the number
% assumptions that the calling code needs to make.
% --
% This class is *NOT* intended to be substituted by a mock class in tests.
% Related class VCDS is however.
%
%
% HOW USING L1R CALIBRATION_TABLE & CALIBRATION_TABLE_INDEX (L1R) WORK
% ====================================================================
% CALIBRATION_TABLE       : CDF L1R global attribute
%   """"Filename of the calibration table(s).""""
%   """"There must as many as entries than the number of calibration table
%   files associated to the L1R file.""""
%
% CALIBRATION_TABLE_INDEX : CDF L1R zVariable
%   """"Index of the calibration table(s) value required to generate L2 data
%   files.""""
%   """"Each CDF record must contain 2 elements: the first element must gives
%   the index of the associated CALIBRATION_TABLE entry (i.e., 0 for the first
%   entry, 1 for the second, etc.). The second element must refer to the index
%   of the value to be used inside the calibration table file.""""
%
% NOTE: Neither exists in L1 datasets.
% NOTE: NBCI is not set (used) for TDS. Therefore no such settings for TDS.
%
% Source: ROC-PRO-DAT-NTT-00006-LES_Iss01_Rev02(ROC_Data_Products).Draft2020-04-06.pdf
%
% Summary
% -------
% L1R GA CALIBRATION_TABLE{CALIBRATION_TABLE_INDEX{iRecord, 1} + 1}
%     == RCT filename
% L1R ZV CALIBRATION_TABLE_INDEX{iRecord, 2}
%     == NBCI
%     == Index/pointer to some calibration value(s) to use in the corresponding
%        RCT. The exact interpretation depends on the RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef VoltageCalibration
  % PROPOSAL: Store all possible versions of TFs internally.
  %   Ex: FTF, ITF, tabulated ITF with extrapolation+interpolation+modification
  %   PRO: Useful for debugging. Can easily inspect & plot FTFs.
  %   NOTE: BIAS & LFR RCTs contain FTFs, TDS RCT contains ITFs.
  %   NOTE: Has to keep track of FTF/ITF before modifications (extrapolation
  %         to 0 Hz, Z=0 for high freq.).
  %   NOTE: Modification (besides inversion) happens on the combined function
  %         which is not stored beforehand.
  %   NOTE: The set of BIAS+LFR+TDS TFs is different from the set of TFs actually
  %         used (combinations of BIAS+LFR and BIAS+TDS respectively).
  %   PROPOSAL: Store LFR TFs as one 1D array of structs with fields: iLsf, iBlts, ftf, itf, ...
  %       PRO: Can easily iterate over.
  %       PRO: For every modification of TFs, can easily add another field for the old
  %            version.
  %       NOTE/CON: All structs/TFs must have the same fields if true struct array.
  %
  %
  %
  % PROPOSAL: Phase out features which have never been used (or not for many
  %           years anyway).
  %   Ex: biasOffsetsDisabled, useBiasTfScalar
  %
  % PROPOSAL: Convert transfer function handles into class(es).
  % PROPOSAL: Move itfHighFreqLimitFraction out of bicas.tf.apply_TF() and
  %           convert it into a modification of the TF instead.
  %   PRO: Simplifies bicas.tf.apply_TF().



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Access=private, Constant)
    % Local TF constants for convenience.
    ONE_TF = @(omegaRps) (ones(size(omegaRps)));
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=public)

    %============================
    % Source of calibration data
    %============================
    Vcds

    %==================================================
    % Settings for what kind of calibration to perform
    %==================================================
    % Correspond to BSO key-values.
    tfMethod
    %
    itfHighFreqLimitFraction
    itfAcConstGainLowFreqRps
    %
    dcDetrendingDegreeOf
    dcRetrendingEnabled
    acDetrendingDegreeOf
    %
    kernelLengthMax
    kernelEdgePolicy
    kernelHannWindowEnabled
    snfEnabled
    snfSubseqMinSamples

    biasOffsetsDisabled

    % What type of BIAS calibration to use.
    useBiasTfScalar



    % Whether to select non-BIAS RCTs using GACT (and NBRI).
    useGactRct

  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Constructor.
    %
    %
    % ARGUMENTS
    % =========
    % Rctdc
    %       Note: Must include BIAS RCTD.
    %
    %
    % NOTES ON INTENDED USAGE
    % =======================
    % The nominal use is that the caller first initializes (argument)
    % RCTDC by
    % (1) loading RCT(s) using
    %     bicas.proc.L1L2.cal.rct.findread.get_nominal_RCTDC(), or
    % (2) manually (for manual debugging/analysis/testing).
    %
    function obj = VoltageCalibration(Vcds, useGactRct, Bso)
      % ASSERTIONS: Arguments
      assert(isa(Vcds, "bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierAbstract"))
      assert(islogical(useGactRct) & isscalar(useGactRct))

      %==============
      % Set obj.Vcds
      %==============
      obj.Vcds = Vcds;

      %======================================================================
      % Store miscellaneous BSO key values
      % ----------------------------------
      % IMPLEMENTATION NOTE: This useful since it:
      %   ** Centralizes the extraction of values from BSO.
      %   ** Is more convenient to access values via shorter field names
      %      (more readable code).
      %   ** Potentially gives faster access to values (better performance).
      %======================================================================
      obj.tfMethod                 = Bso.get_fv('PROCESSING.CALIBRATION.TF.METHOD');

      obj.itfHighFreqLimitFraction = Bso.get_fv('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION');
      % NOTE: Converts Hz-->rad/s
      obj.itfAcConstGainLowFreqRps = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_CONST_GAIN_LOW_FREQ_HZ') * 2*pi;

      obj.dcDetrendingDegreeOf     = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE');
      obj.dcRetrendingEnabled      = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED');
      obj.acDetrendingDegreeOf     = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE');

      obj.kernelLengthMax          = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.MAX_LENGTH_SAMPLES');
      obj.kernelEdgePolicy         = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.EDGE_POLICY');
      obj.kernelHannWindowEnabled  = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.HANN_WINDOW_ENABLED');
      obj.snfEnabled               = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED');
      obj.snfSubseqMinSamples      = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES');

      obj.biasOffsetsDisabled      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED');

      %-------------------------
      % Set obj.useBiasTfScalar
      %-------------------------
      settingBiasTf                = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
      switch(settingBiasTf)
        case 'FULL'
          obj.useBiasTfScalar = false;
        case 'SCALAR'
          obj.useBiasTfScalar = true;
        otherwise
          error(...
            'BICAS:Assertion:ConfigurationBug', ...
            ['Illegal value for setting',...
            ' PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF="%s".'], ...
            settingBiasTf)
      end


      %============================
      % Store some argument values
      %============================
      obj.useGactRct = useGactRct;
    end



    % Calibrate all voltages. Function will choose the more specific algorithm
    % internally.
    %
    %
    % ARGUMENTS
    % =========
    % NbriFpa, NbciFpa
    %       NOTE: Only one record of NBRI, NBCI!.
    %
    function bltsSamplesAvoltCa = calibrate_voltage_TM_to_avolt(obj, ...
        dtSec, bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, NbriFpa, NbciFpa)

      % ASSERTIONS
      assert(iscell(bltsSamplesTmCa))
      irf.assert.sizes(...
        dtSec,           [-1], ...
        bltsSamplesTmCa, [-1])
      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))
      assert(islogical(isTdsCwf) && isscalar(isTdsCwf))



      %===============================================================
      % Obtain detrending/retrending settings for bicas.tf.apply_TF()
      %===============================================================
      if bicas.proc.L1L2.const.SSID_is_AC(CalSettings.ssid)
        % IMPLEMENTATION NOTE: DC is (optionally) detrended via
        % bicas.tf.apply_TF() in the sense of a fit being removed, TF applied,
        % and the fit then added back. That same algorithm, or at least adding
        % back the fit, is by its nature inappropriate for non-lowpass filters,
        % i.e. for AC. (The fit can not be scaled with the 0 Hz signal
        % amplitude)

        % NOTE: Using setting for AC specifically!
        detrendingDegreeOf = obj.acDetrendingDegreeOf;
        retrendingEnabled  = false;                   % NOTE: HARDCODED SETTING.
      else
        % NOTE: Using settings for DC specifically!
        detrendingDegreeOf = obj.dcDetrendingDegreeOf;
        retrendingEnabled  = obj.dcRetrendingEnabled;
      end



      %=============================
      % Obtain itfAvpt, offsetAvolt
      %=============================
      [itfAvpt, offsetAvolt] = obj.get_voltage_calibration_ITF_offset(...
        isLfr       = isLfr, ...
        isTdsCwf    = isTdsCwf, ...
        CalSettings = CalSettings, ...
        NbriFpa     = NbriFpa, ...
        NbciFpa     = NbciFpa);



      %======================================
      % Apply BIAS+LFR/TDS+offset to samples
      %======================================
      bltsSamplesAvoltCa = cell(size(bltsSamplesTmCa));    % Pre-allocate
      for i = 1:numel(bltsSamplesTmCa)
        tempSamplesAvolt = bicas.tf.apply_TF(...
          dtSec(i), ...
          bltsSamplesTmCa{i}, ...
          itfAvpt, ...
          method                  = obj.tfMethod, ...
          detrendingDegreeOf      = detrendingDegreeOf, ...
          retrendingEnabled       = retrendingEnabled, ...
          tfHighFreqLimitFraction = obj.itfHighFreqLimitFraction, ...
          kernelLengthMax         = obj.kernelLengthMax, ...
          kernelEdgePolicy        = obj.kernelEdgePolicy, ...
          kernelHannWindowEnabled = obj.kernelHannWindowEnabled, ...
          snfEnabled              = obj.snfEnabled, ...
          snfSubseqMinSamples     = obj.snfSubseqMinSamples);

        bltsSamplesAvoltCa{i, 1} = tempSamplesAvolt + offsetAvolt;
      end

    end



    % Function for producing the exact ITF used by BICAS for all calibration.
    %
    % RATIONALE FOR EXISTENCE
    % =======================
    % This function exists to make it possible to manually reproduce and
    % inspect the exact TF+offset used by BICAS, from outside of BICAS proper.
    % NOTE: This still excludes information on detrending-retrending, use of
    % Hann window, kernel length.
    %
    function [itfAvpt, offsetAvolt] = get_voltage_calibration_ITF_offset(obj, A)
      % PROPOSAL: Redefine to return all information (arguments) used for
      %           calling bicas.tf.apply_TF(), including detrending-retrending,
      %           use of Hann window, kernel length, SNF options.
      %   PRO: Potentially useful for the same reasons as in the rationale for
      %        this function.
      %     CON: Not sure if it is important yet. Issues associated with other
      %          "calibration parameters" are less subtle than for TFs/kernels.
      %   CON: Many return values. Need to use struct/class.
      %
      % PROPOSAL: Reorg. code to have one function for returning the modified
      % TF as now done in bicas.tf.apply_TF(), i.e. remove high frequencies.

      arguments
        obj
        A.isLfr
        A.isTdsCwf
        A.CalSettings
        A.NbriFpa
        A.NbciFpa
      end

      assert(islogical(A.isTdsCwf) && isscalar(A.isTdsCwf))
      assert(isa(A.NbriFpa, "bicas.utils.FPArray"))
      assert(isa(A.NbciFpa, "bicas.utils.FPArray"))

      iBlts = A.CalSettings.iBlts;



      % Normalize NbriFpa.
      if ~obj.useGactRct
        % Emulating having an NBRI value.
        A.NbriFpa = bicas.utils.FPArray(1, 'NO_FILL_POSITIONS');
      end



      if bicas.proc.L1L2.const.SSID_is_ASR(A.CalSettings.ssid)
        % CASE: ASR SSID

        %========================
        % Obtain BIAS ITF+offset
        %========================
        [itfAvpiv, kItfAvpiv, offsetAvolt] = obj.Vcds.get_BIAS_ITF_and_offset(...
          A.CalSettings.ssid, ...
          A.CalSettings.isAchg, ...
          A.CalSettings.iCalibTimeL, ...
          A.CalSettings.iCalibTimeH);

        if obj.useBiasTfScalar
          % NOTE: Overwrites pre-existing "itfAvpiv" which still needs to be
          % created first in order to obtain kItfAvpiv (needed here).
          itfAvpiv = @(omegaRps) (ones(size(omegaRps)) * kItfAvpiv);
        end

        %====================
        % Obtain LFR/TDS ITF
        %====================
        if A.isLfr
          %===========
          % CASE: LFR
          %===========
          itfIvpt = obj.Vcds.get_LFR_ITF(iBlts, A.NbriFpa, A.NbciFpa, A.CalSettings.iLsf);
        else
          %===========
          % CASE: TDS
          %===========
          if A.isTdsCwf
            % CASE: TDS CWF
            itfIvpt = obj.Vcds.get_TDS_CWF_ITF(iBlts, A.NbriFpa, A.NbciFpa);
          else
            % CASE: TDS RSWF
            itfIvpt = obj.Vcds.get_TDS_RSWF_ITF(iBlts, A.NbriFpa, A.NbciFpa);
          end
        end

        itfAvpt = obj.combine_BIAS_ITF_and_LFR_TDS_ITF(...
          itfIvpt, itfAvpiv, ...
          bicas.proc.L1L2.const.SSID_is_AC(A.CalSettings.ssid), ...
          obj.itfAcConstGainLowFreqRps);

      else
        % CASE: Non-ASR SSID
        itfAvpt     = obj.ONE_TF;
        offsetAvolt = 0;
      end

      if obj.biasOffsetsDisabled
        % NOTE: Overwrites pre-existing value of "offsetAvolt".
        offsetAvolt = 0;
      end
    end



    function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
      iCalibL = obj.Vcds.get_BIAS_calibration_time_index_L(Epoch);
    end



    function iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)
      iCalibH = obj.Vcds.get_BIAS_calibration_time_index_H(Epoch);
    end



  end    % methods(Access=private)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static, Access=private)



    % Combine BIAS and LFR/TDS ITFs, with a special modification for (LFR) AC.
    %
    function itf = combine_BIAS_ITF_and_LFR_TDS_ITF(...
        itfLfr, itfBias, isAc, acConstGainLowFreqRps)

      assert(isscalar(isAc), islogical(isAc))

      itf = @(omegaRps) (TF_product(omegaRps));

      if isAc
        % NOTE: Modifies the already created, combined LFR+BIAS TF.

        zLimit = itf(acConstGainLowFreqRps);

        itf = @(omegaRps) (bicas.proc.L1L2.cal.utils.TF_LF_constant_abs_Z(...
          itf, omegaRps, acConstGainLowFreqRps, zLimit));
      end


      %###################################################################
      % IMPLEMENTATION NOTE: In principle, this function is quite
      % unnecessary for multiplying TFs, but it is useful for putting
      % breakpoints in when debugging TFs which are built from layers of
      % anonymous functions and function handles.
      function Z = TF_product(omegaRps)
        Z = itfLfr(omegaRps) ...
          .* ...
          itfBias(omegaRps);
      end
    end



  end    % methods(Static, Access=private)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static, Access=public)



    % TODO-NI: Where does the parasitic capacitance TF fit into the calibration formulas?
    % TODO-NI: What parasitic capacitance value(s) should one use?
    % PROPOSAL: Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
    %
    %         function tfZ = parasitic_capacitance_TF(tfOmega)
    %             % Calculate Z(omega) values for TF representing parasitic
    %             % capacitances (based on analytic function).
    %
    %             % Function name? "Input capacitance"?
    %             % Not read R & C from constants here? Submit as arguments?
    %             capacitanceFarad =
    %             impedanceOhm     =
    %
    %             % Correct for a TF?!
    %             % 2020-11-11, EJ: Not same as in paper note.
    %             tfZ = 1 / (1 + 1j*tfOmega*capacitanceFarad*impedanceOhm);
    %
    %             error('BICAS:OperationNotImplemented', 'Function not implemented Yet.')
    %         end



  end    % methods(Static, Access=public)



end
