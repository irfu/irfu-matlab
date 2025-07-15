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
% --
% BUG: Can likely not handle data with non-ASR SSID (Unknown, GND, or 2.5V
%      Ref).
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
% NOTE: ZVCTI2 is not set (used) for TDS. Therefore no such settings for TDS.
%
% Source: ROC-PRO-DAT-NTT-00006-LES_Iss01_Rev02(ROC_Data_Products).Draft2020-04-06.pdf
%
% Summary
% -------
% L1R GA CALIBRATION_TABLE{CALIBRATION_TABLE_INDEX{iRecord, 1} + 1}
%     == RCT filename
% L1R ZV CALIBRATION_TABLE_INDEX{iRecord, 2}
%     == ZVCTI2
%     == Index/pointer to some calibration value(s) to use in the corresponding
%        RCT. The exact interpretation depends on the RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef VoltageCalibration
  % PROPOSAL: Separate subclasses for different types of data. At least separate
  %           for LFR, TDS-CWF, TDS-RSWF.
  %   PROPOSAL: Methods
  %     bltsSamplesAVoltCa = get_LFR_ITF(obj, ...
  %     bltsSamplesAVoltCa = get_TDS_CWF_ITF(obj, ...
  %     bltsSamplesAVoltCa = get_TDS_RSWF_ITF(obj, ...
  %   should then be separate implementations of one abstract superclass mathod.
  %     PROPOSAL: Separate non-abstract superclass method (which handles
  %       allVoltageCalibDisabled, ufv, useGact(=false)) calls the subclass
  %        method.
  %   PROPOSAL: Separate subclass for ignoring calibration.
  %     CON: There is too much code that is not nominal ITFs (from RCTs), that
  %          could need testing and should not be "mocked away".
  %   PROPOSAL: Separate subclass for mocking in automated tests.
  %   CON/PROBLEM: Methods duplicate the derivation of BIAS TF (by calling one
  %                function) and the combining of BIAS TF & LFR/TDS TF.
  %                Combining TFs should really be done centrally, and obtaining
  %                BIAS TF should be done once. ==> Separate methods/sources
  %                for producing BIAS TF and LFR/TDS TF. ==> Ambiguous how to
  %                apply "mocking design pattern". Does not make sense for
  %                subclass to override both or only one of these methods.
  %     PROPOSAL: Test subclass overrides both BIAS TF and LFR/TDS TF methods.
  %       PRO: Simple.
  %         PRO: Realistically, a test always needs to override both BIAS TF and
  %              LFR/TDS TFs at the same time.
  %       CON: Conceptually ugly.
  %       CON/PROBLEM: Test subclass is "forced" to load RCTDs!!!!!
  %         NOTE: The BIAS/LFR/TDS TF-producing methods need RCTDs but the
  %               implementing subclasses can load these.
  %         PRO: get_BIAS_calibration_time_index_L/H() need BIAS RCTDs.
  %              ==> Re-implemented in multiple subclasses, unless the
  %              superclass (needed for testing) should also load the BIAS RCT!
  %           NOTE: CCAL also needs BIAS RCTD for implementing
  %                 get_BIAS_calibration_time_index_L().
  %       PRO: Can then abolish bicas.proc.L1L2.cal.VoltageCalibrationAbstract. -- IMPLEMENTED
  %     --
  %     PROPOSAL: Class uses separately mockable classes for producing BIAS TF
  %               and LFR/TDS TF.
  %       CON: Gets too many classes
  %         VCAL (no abstract+impl.+test)
  %           Non-BIAS TF: Abstract+test + LFR, TDS-CWF, TDS-RSWF (impl.)
  %           BIAS TF:     Abstract+test + impl.
  %         CON: Can (probably) abolish VCAL abstract+test.
  %         CON-PROPOSAL: Abolish VCAL?!!!
  %           PRO: The only "real" VCAL function is calibrate_voltage_TM_to_avolt().
  %       CON: Conceptually, LFR+TDS calibration is not exchangeable (can not be
  %            substituted for each other) and is applied to different types of
  %            data (with different flags).
  %            ==> Should not have separate subclasses.
  %       PROPOSAL: One abstract+impl+test set of classes for BIAS+LFR/TDS TFs. -- IMPLEMENTED
  %       PROPOSAL: Merge CCAL+BIAS TF producer.
  %         PRO: Both use BIAS RCT.
  %       PROPOSAL: Class names
  %         Offset
  %         ~TF
  %           CON: Not all calibration uses a (frequency-dependent) TF.
  %         Calibration
  %         Voltage
  %         BIAS, LF/TDS
  %         --
  %         BIAS
  %           BCAL=BiasCalibration
  %           BiasUnitCalibration
  %           VCBH=VoltageCalibrationBiasHw
  %         NonBias
  %           VoltageCalibrationNonBias
  %           VCLT=VoltageCalibrationLfrTdsAbstract
  %             VoltageCalibrationLfr
  %             VoltageCalibrationTdsCwf
  %             VoltageCalibrationTdsRswf
  %           TmBltsVoltageCalibration
  %
  %       NOTE: BIAS+LFR/TDS TF classes must separately load RCTDs.
  %         PRO?: Handling of BIAS and LFR/TDS RCTs are different anyway.
  %     --
  %     PROPOSAL: VCAL makes it really easy for one subclass to combine LFR/TDS
  %               TF with BIAS TF by supplying method for it.
  %
  %
  %
  % PROPOSAL: Tests for non-ASR SSIDs.
  %
  %
  %
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
  % PROPOSAL: Phase out some features which could potentially be achieved by
  %           using alternative subclasses.
  %   Ex: allVoltageCalibDisabled
  %   NOTE: May want to disable calibration (use unit TF) for non-ASR signals.
  %     NOTE: Disable reconstruction for this case too?!!
  %   PROPOSAL: Convert allVoltageCalibDisabled to method argument (currently
  %             set at initialization).
  % PROPOSAL: Phase out features which have never been used (or not for many
  %           years anyway).
  %   Ex: lfrTdsTfDisabled
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

    %==================
    % Calibration data
    %==================
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
    kernelEdgePolicy
    kernelHannWindow
    snfEnabled
    snfSubseqMinSamples

    % What type of calibration to use.
    allVoltageCalibDisabled    % Use TM values (not set to NaN).
    useBiasTfScalar
    lfrTdsTfDisabled

    % Whether to select non-BIAS RCTs using GACT (and ZVCTI).
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
    % RCTDC
    % (1) by loading RCTs using
    %     bicas.proc.L1L2.cal.rct.findread.find_read_nonBIAS_RCTs_by_regexp(),
    % (2) by loading RCT(s) using
    %     bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_BRVF_and_ZVCTI_GACT()
    % or
    % (3) manually (for manual debugging/analysis/testing).
    %
    function obj = VoltageCalibration(Vcds, useGactRct, Bso)
      % ASSERTIONS: Arguments
      assert(isa(Vcds, "bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierAbstract"))
      assert(islogical(useGactRct) & isscalar(useGactRct))

      %==============
      % Set obj.Vcds
      %==============
      obj.Vcds = Vcds;

      %===============================================================
      % Store miscellaneous BSO key values
      % ----------------------------------
      % IMPLEMENTATION NOTE: This useful since it is:
      %   ** More convenient to access values via shorter field names
      %      (more readable code).
      %   ** Potentially gives faster access to values (better
      %      performance).
      %===============================================================
      obj.tfMethod                 = Bso.get_fv('PROCESSING.CALIBRATION.TF.METHOD');

      obj.itfHighFreqLimitFraction = Bso.get_fv('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION');
      % NOTE: Converts Hz-->rad/s
      obj.itfAcConstGainLowFreqRps = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_CONST_GAIN_LOW_FREQ_HZ') * 2*pi;

      obj.dcDetrendingDegreeOf     = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE');
      obj.dcRetrendingEnabled      = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED');
      obj.acDetrendingDegreeOf     = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE');

      obj.kernelEdgePolicy         = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.EDGE_POLICY');
      obj.kernelHannWindow         = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.HANN_WINDOW_ENABLED');
      obj.snfEnabled               = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED');
      obj.snfSubseqMinSamples      = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES');

      obj.allVoltageCalibDisabled  = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLED');
      obj.lfrTdsTfDisabled         = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED');

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
    % zvcti
    %       NOTE: Only one record (row) of ZVCTI! Not entire ZV.
    % ufv
    %       Scalar logical.
    %       Whether to set output voltages to NaN (representing fill values)
    %       and thus not execute any (real) calibration.
    %       RATIONALE: This option is useful to
    %       (1) potentially speed up BICAS when it is known that
    %           data will be overwritten with fill values later.
    %       (2) avoid executing calibration algorithms when it is
    %           known that there is no calibration configuration anyway
    %           Ex: LFR zVar BW=0 ==> zvcti(1,:) value is illegal.
    %               ==> Can not calibrate.
    %           Note: This means that this function technically accepts
    %           an illegal calibration configuration when this argument is set
    %           to true.
    %
    function bltsSamplesAVoltCa = calibrate_voltage_TM_to_avolt(obj, ...
        dtSec, bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, zvcti, ufv)

      % ASSERTIONS
      assert(iscell(bltsSamplesTmCa))
      irf.assert.sizes(...
        zvcti,           [ 1, 2], ...
        dtSec,           [-1], ...
        bltsSamplesTmCa, [-1])
      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))
      assert(islogical(ufv)      && isscalar(ufv))
      assert(islogical(isTdsCwf) && isscalar(isTdsCwf))

      iBlts = CalSettings.iBlts;
      iLsf  = CalSettings.iLsf;



      % Set iNonBiasRct.
      if obj.useGactRct
        % NOTE: Incrementing by one (index into MATLAB array).
        iNonBiasRct = 1 + zvcti(1, 1);
      else
        % Emulating having a ZVCTI value.
        iNonBiasRct = 1;
      end
      % NOTE: NOT incrementing value by one, since the variable's meaning
      % can vary between LFR, TDS-CWF, TDS-RSWF.
      zvcti2 = zvcti(1,2);



      %=========================================
      % Obtain settings for bicas.tf.apply_TF()
      %=========================================
      if bicas.proc.L1L2.const.SSID_is_AC(CalSettings.ssid)
        % IMPLEMENTATION NOTE: DC is (optionally) detrended via
        % bicas.tf.apply_TF() in the sense of a linear fit being removed, TF
        % applied, and then added back. That same algorithm, or at least adding
        % back the fit, is by its nature inappropriate for non-lowpass filters,
        % i.e. for AC. (The fit can not be scaled with the 0 Hz signal
        % amplitude)
        detrendingDegreeOf = obj.acDetrendingDegreeOf;    % NOTE: AC!
        retrendingEnabled  = false;                   % NOTE: HARDCODED SETTING.
      else
        detrendingDegreeOf = obj.dcDetrendingDegreeOf;    % NOTE: DC!
        retrendingEnabled  = obj.dcRetrendingEnabled;
      end



      %================
      % Obtain itfAvpt
      %================
      if obj.allVoltageCalibDisabled || ufv

        if obj.allVoltageCalibDisabled
          % CASE: Set voltages to TM values.

          itfAvpt     = obj.ONE_TF;
          offsetAvolt = 0;
        end
        if ufv
          % CASE: Set voltages to NaN.

          % IMPLEMENTATION NOTE: Potentially overwrites value set in above "if"
          % statement.
          itfAvpt     = bicas.const.NAN_TF;
          offsetAvolt = NaN;
        end

      else
        %========================
        % Obtain BIAS ITF+offset
        %========================
        [itfAvpiv, offsetAvolt] = obj.Vcds.get_BIAS_ITF_and_offset(...
          CalSettings.ssid, ...
          CalSettings.isAchg, ...
          CalSettings.iCalibTimeL, ...
          CalSettings.iCalibTimeH);

        %====================
        % Obtain LFR/TDS ITF
        %====================
        if obj.lfrTdsTfDisabled
          itfIvpt = obj.ONE_TF;
        else
          if isLfr
            %===========
            % CASE: LFR
            %===========
            itfIvpt = obj.Vcds.get_LFR_ITF(iBlts, iLsf, iNonBiasRct, zvcti2);
          else
            %===========
            % CASE: TDS
            %===========
            if isTdsCwf
              % CASE: TDS CWF
              itfIvpt = obj.Vcds.get_TDS_CWF_ITF(iBlts, iNonBiasRct, zvcti2);
            else
              % CASE: TDS RSWF
              itfIvpt = obj.Vcds.get_TDS_RSWF_ITF(iBlts, iNonBiasRct, zvcti2);
            end
          end
        end

        itfAvpt = bicas.proc.L1L2.cal.utils.combine_BIAS_ITF_and_LFR_TDS_ITF(...
          itfIvpt, itfAvpiv, ...
          bicas.proc.L1L2.const.SSID_is_AC(CalSettings.ssid), ...
          obj.itfAcConstGainLowFreqRps);
      end



      %======================================
      % Apply BIAS+LFR/TDS+offset to samples
      %======================================
      bltsSamplesAVoltCa = cell(size(bltsSamplesTmCa));
      for i = 1:numel(bltsSamplesTmCa)
         tempSamplesAVolt = bicas.tf.apply_TF(...
          dtSec(i), ...
          bltsSamplesTmCa{i}(:), ...
          itfAvpt, ...
          'method',                  obj.tfMethod, ...
          'detrendingDegreeOf',      detrendingDegreeOf, ...
          'retrendingEnabled',       retrendingEnabled, ...
          'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction, ...
          'kernelEdgePolicy',        obj.kernelEdgePolicy, ...
          'kernelHannWindow',        obj.kernelHannWindow, ...
          'snfEnabled',              obj.snfEnabled, ...
          'snfSubseqMinSamples',     obj.snfSubseqMinSamples);

         bltsSamplesAVoltCa{i, 1} = tempSamplesAVolt + offsetAvolt;
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
