%
% Nominal implementation of superclass for calibrating voltages.
%
% An instance of this class contains
%   (1) relevant settings (loaded from BSO) on how to calibrate voltages, and
%   (2) actual calibration data (from RCTs).
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
%
% Source: ROC-PRO-DAT-NTT-00006-LES_Iss01_Rev02(ROC_Data_Products).Draft2020-04-06.pdf
%
% Summary
% -------
% GA CALIBRATION_TABLE{CALIBRATION_TABLE_INDEX{iRecord, 1} + 1}
%     == RCT filename
% ZV CALIBRATION_TABLE_INDEX{iRecord, 2}
%     == ZVCTI2
%     == Index/pointer to some calibration value(s) to use in the corresponding
%        RCT. The exact interpretation depends on the RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef VoltageCalibrationImpl < bicas.proc.L1L2.cal.VoltageCalibrationAbstract
  % All methods as of 2025-07-14
  % ----------------------------
  % function obj = VoltageCalibrationImpl(Rctdc, useGactRct, useZvcti2, Bso)
  % function bltsSamplesAVoltCa = calibrate_voltage_all(obj, ...
  % function [itfAvpiv, offsetAvolt] = get_BIAS_ITF(obj, ...
  % function itfIvpt = get_LFR_ITF(obj, iBlts, iLsf, iNonBiasRct, zvcti2)
  % function itfIvpt = get_TDS_CWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)
  % function itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)
  % --
  % function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
  % function iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)
  %
  %
  %
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
  %       CON/PROBLEM: Test subclass is forced to load RCTDs!!!!!
  %         NOTE: The BIAS/LFR/TDS TF-producing methods need RCTDs but the
  %               implementing subclasses can load these.
  %         PRO: get_BIAS_calibration_time_index_L/H() need BIAS RCTDs.
  %           NOTE: CCAL also needs BIAS RCTD for implementing
  %                 get_BIAS_calibration_time_index_L().
  %       NOTE: Can then abolish bicas.proc.L1L2.cal.VoltageCalibrationAbstract.
  %     --
  %     PROPOSAL: Class uses separately mockable classes for producing BIAS TF
  %               and LFR/TDS TF.
  %       CON: Gets too many classes
  %         VCAL (no abstract+impl.+test)
  %           Non-BIAS TF: Abstract+test + LFR, TDS-CWF, TDS-RSWF (impl.)
  %           BIAS TF:     Abstract+test + impl.
  %         CON: Can (probably) abolish VCAL abstract+test.
  %         CON-PROPOSAL: Abolish VCAL?!!!
  %           PRO: The only "real" VCAL function is calibrate_voltage_all().
  %       PROPOSAL: Merge CCAL+BIAS TF producer.
  %         PRO: Both use BIAS RCT.
  %         PROPOSAL: Name
  %           BCAL=BiasCalibration, BiasUnitCalibratio, BiasHwCalibration.
  %       NOTE: BIAS+LFR/TDS TF classes must separately load RCTDs.
  %         PRO?: Handling of BIAS and LFR/TDS RCTs are different anyway.
  %     --
  %     PROPOSAL: VCAL makes it really easy for one subclass to combine LFR/TDS
  %               TF with BIAS TF.
  %   --
  %   PRO: Large class: 1004 rows. /2025-07-14
  %
  %
  %
  % PROPOSAL: Refactor to facilitate automated testing.
  %   PROBLEM: Though it does not read RCTs, the corresponding data
  %            structures are complex and would be hard to create test data for(?)
  %   NOTE: Uses BSO, and many of its values.
  %   PROBLEM: Calls complex function that should (primarily) be tested separately:
  %            bicas.tf.apply_TF_freq()
  %       PROPOSAL: Replace with function handle, set in constructor from
  %                 argument ("dependency injection"). Unit tests can then mock
  %                 bicas.tf.apply_TF_freq().
  %           CON?: Introduces more interface (arguments) only due to testing.
  %       PROPOSAL: Do automatic testing by having the tests call
  %                 bicas.tf.apply_TF_freq() to generate results to compare with.
  %           CON: Relies on the implementation of what is being tested.
  %           CON: Can not test arguments sent to bicas.tf.apply_TF_freq().
  %               CON: Relies on the implementation of what is being tested.
  %       NOTE: (1) As a function/code module,
  %                   bicas.proc.L1L2.cal.VoltageCalibrationImpl
  %                   encloses/contains/"secretly uses"
  %                   bicas.tf.apply_TF_freq().
  %             (2) For testing, one wants to verify the path (both ways) between
  %                   bicas.proc.L1L2.cal.VoltageCalibrationImpl and
  %                   bicas.tf.apply_TF_freq().
  %             ==> One wants to test one unit of code at a time, but what a
  %                "unit" is ambiguous:
  %                 one wants small units of code
  %                 unit is ambiguous when a unit uses/call other unit(s).
  %
  %
  %
  % TODO-DEC: How distribute the calibration formulas/algorithms between
  %   (1) calibrate_* functions,
  %   (2) functions that select relevant calibration data
  %   (2) RCT reading functions,
  %   (3) classes that store TFs?
  %   --
  %   Ex: Invert the (rat.func., tabulated) TFs
  %   Ex: Extrapolate the tabulated TFs to zero
  %   Ex: Extrapolate the tabulated LFR TF to higher frequencies.
  %   Ex: Modify AC TF at lower frequencies (change amplitude, keep phase)
  %   Ex: Interpolation algorithm of tabulated TFs.
  %   Ex: If one modifies the TFs before applying the inverse (hypothetical; not implemented)
  %   --
  %   NEED: Plot all TFs.
  %   NEED: Plot all versions of any particular TF as it gets modified.
  %   NEED: Plot all TFs used for a particular calibration case
  %         (when calibrating using bicas.caib only, without BICAS).
  %   PROPOSAL: Separate modifications/choices/code that
  %       (1) only need to be done/run once during the execution:
  %               modification of calibration data,
  %       (2) are done every calibration (call to calibrate_*):
  %               exact algorithms/formulas through which data is fed
  %   PROPOSAL: read_*_RCT should not modify any calibration data, just store it:
  %             Not invert TFs, not extrapolate TFs to 0 Hz.
  %       CON: Must separate unmodified and modified calibration data.
  %       CON: "Must" label variables that they are modified/unmodified.
  %       CON-PROBLEM: No clear line for what is a modification or not.
  %           NOTE: Difference between
  %               (1) modification (information potentially destroyed), and
  %               (2) conversion (no information destroyed; e.g. format conversion)
  %           Ex: TF frequency Hz-->omega rad/s
  %           Ex: TF amplitude+phase-->Z
  %           Ex: Apply upper frequency cut-off to ITF, in particular analytical ITFs.
  %           Ex: Extrapolate tabulated TF
  %           Ex: How/where make different choices for how to calibrate?
  %               (1) No calibration
  %               (2) Scalar calibration (a) with/(b) without offsets
  %               (3) Full calibration
  %               (4) Full calibration except without parasitic capacitance.
  %       TODO-DEC: Where is it natural to modify calibration data then?!
  %   PROPOSAL: General philosophy should be that calibrate_* chooses as much as
  %             possible, and thus chooses different functions to call.
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
  % BUG: Can likely not handle data with SSID = Unknown or 2.5V Ref, at least
  %      not for LFR.
  %   PROPOSAL: Tests.
  %
  %
  %
  % PROPOSAL: Phase out some features which could potentially be achieved by
  %           using alternative subclasses.
  %   Ex: ufv, allVoltageCalibDisabled
  % PROPOSAL: Phase out features which have never been used (or not for many
  %           years anyway).
  %   Ex: lfrTdsTfDisabled
  %
  % PROPOSAL: calibrate_voltage_all() should centralize
  %   (1) one shared loop over sample sequences in argument cell array,  -- IMPLEMENTED
  %   (2) combining BIAS TF and LFR/TDS TF,                              -- IMPLEMENTED
  %   (3) call to bicas.tf.apply_TF(),                                   -- IMPLEMENTED
  %   (4) add BIAS offset                                                -- IMPLEMENTED
  %   rather than having the type-specific methods implement it separately:
  %   get_LFR_ITF/TDS_CWF/TDS_RSWF().
  %   Make those function return a transfer function instead.
  %
  %
  %
  % PROPOSAL: Convert transfer function handles into class(es).



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Access=private, Constant)
    % Local TF constants for convenience.
    NAN_TF = @(omegaRps) (omegaRps * NaN);
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

    % RCT calibration data
    Rctdc;

    % Non-RCT calibration data
    % ------------------------
    % BIAS scalar (simplified) calibration, not in the RCTs. For
    % debugging/testing purposes.
    BiasScalarGain



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
    biasOffsetsDisabled
    lfrTdsTfDisabled

    % Whether to select non-BIAS RCTs using GACT (and ZVCTI).
    useGactRct
    % Whether to use ZVCTI2 for calibration.
    useZvcti2

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
    %
    % IMPLEMENTATION NOTE
    % ===================
    % The class (instance methods, including constructor) deliberately does
    % not itself read the RCTs, nor figure out which ones should be read.
    % This is useful since
    % ** it completely separates
    %       (a) algorithms for determining RCTs to load, and
    %       (b) reading RCT,
    %    from the class (better modularization, better for automatic test
    %    code).
    % ** it makes it possible to inspect & modify the RCT content before
    %    submitting it to bicas.proc.L1L2.cal.VoltageCalibrationAbstract
    % ** it simplifies the constructor.
    %
    function obj = VoltageCalibrationImpl(Rctdc, useGactRct, useZvcti2, Bso)

      % ASSERTIONS: Arguments
      assert(islogical(useGactRct) & isscalar(useGactRct))
      assert(islogical(useZvcti2)  & isscalar(useZvcti2))
      % Rctdc
      assert(isa(Rctdc, 'bicas.proc.L1L2.cal.RctdCollection'))



      %===============
      % Set obj.Rctdc
      %===============
      obj.Rctdc = Rctdc;



      %==================================================================
      % Store miscellaneous BSO key values
      % ----------------------------------
      % IMPLEMENTATION NOTE: This useful since it is:
      %   ** More convenient to access values via shorter field names
      %      (more readable code).
      %   ** Potentially gives faster access to values (better
      %      performance).
      %==================================================================
      obj.BiasScalarGain.alphaIvpav      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
      obj.BiasScalarGain.betaIvpav       = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
      obj.BiasScalarGain.gammaIvpav.achg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
      obj.BiasScalarGain.gammaIvpav.aclg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');

      obj.tfMethod                       = Bso.get_fv('PROCESSING.CALIBRATION.TF.METHOD');

      obj.itfHighFreqLimitFraction       = Bso.get_fv('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION');
      % NOTE: Converts Hz-->rad/s
      obj.itfAcConstGainLowFreqRps       = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_CONST_GAIN_LOW_FREQ_HZ') * 2*pi;

      obj.dcDetrendingDegreeOf           = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE');
      obj.dcRetrendingEnabled            = Bso.get_fv('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED');
      obj.acDetrendingDegreeOf           = Bso.get_fv('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE');

      obj.kernelEdgePolicy               = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.EDGE_POLICY');
      obj.kernelHannWindow               = Bso.get_fv('PROCESSING.CALIBRATION.TF.KERNEL.HANN_WINDOW_ENABLED');
      obj.snfEnabled                     = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED');
      obj.snfSubseqMinSamples            = Bso.get_fv('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES');

      obj.allVoltageCalibDisabled        = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLED');
      obj.biasOffsetsDisabled            = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED');
      obj.lfrTdsTfDisabled               = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED');

      %-------------------------
      % Set obj.useBiasTfScalar
      %-------------------------
      settingBiasTf                      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
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
      obj.useZvcti2  = useZvcti2;
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
    function bltsSamplesAVoltCa = calibrate_voltage_all(obj, ...
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
          itfAvpt     = obj.NAN_TF;
          offsetAvolt = NaN;
        end

      else
        %========================
        % Obtain BIAS ITF+offset
        %========================
        [itfAvpiv, offsetAvolt] = obj.get_BIAS_ITF(...
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
            itfIvpt = obj.get_LFR_ITF(iBlts, iLsf, iNonBiasRct, zvcti2);
          else
            %===========
            % CASE: TDS
            %===========
            if isTdsCwf
              % CASE: TDS CWF
              itfIvpt = obj.get_TDS_CWF_ITF(iBlts, iNonBiasRct, zvcti2);
            else
              % CASE: TDS RSWF
              itfIvpt = obj.get_TDS_RSWF_ITF(iBlts, iNonBiasRct, zvcti2);
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



    % Return BIAS ITF and BIAS offset:
    %
    % NOTE: May return calibration values corresponding to scalar
    % calibration, depending on BSO:
    %
    %
    % ARGUMENTS
    % =========
    % isAchg
    %       NUMERIC value: 0=Off, 1=ON, or NaN=Value not known.
    %       IMPLEMENTATION NOTE: Needs value to represent that isAchg
    %       is unknown. Sometimes, if isAchg is unknown, then it is
    %       useful to process as usual since some of the data can still be
    %       derived/calibrated, so that the caller does not need to handle
    %       the special case.
    %
    function [itfAvpiv, offsetAvolt] = get_BIAS_ITF(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)

      % PROPOSAL: Log warning message when simultaneously isAchg=NaN
      %           and the value is needed.
      % PROPOSAL: Separate functions for TF and scalar calibration.

      % ASSERTION
      assert(bicas.proc.L1L2.const.is_SSID(ssid) & isscalar(ssid))
      assert(bicas.proc.L1L2.const.SSID_is_ASR(ssid))
      assert(isscalar(isAchg) && isnumeric(isAchg))
      assert(isscalar(iCalibTimeL))
      assert(isscalar(iCalibTimeH))

      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');
      BiasRctd   = BiasRctdCa{1};

      %###################################################################
      % kFtfIvpav = Multiplication factor "k" that represents/replaces the
      % (forward) transfer function for scalar calibration.
      asid         = bicas.proc.L1L2.const.SSID_ASR_to_ASID( ssid);
      asidCategory = bicas.proc.L1L2.const.get_ASID_category(asid);
      antennas     = bicas.proc.L1L2.const.get_ASID_antennas(asid);
      switch(asidCategory)

        case 'DC_SINGLE'

          % NOTE: List of ITFs for different times.
          itfAvpiv    = BiasRctd.ItfSet.dcSingleAvpiv{iCalibTimeL};
          kFtfIvpav   = obj.BiasScalarGain.alphaIvpav;
          offsetAvolt = BiasRctd.dcSingleOffsetsAVolt(iCalibTimeH, antennas);

        case 'DC_DIFF'

          itfAvpiv = BiasRctd.ItfSet.dcDiffAvpiv{iCalibTimeL};
          kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
          if     isequal(antennas, [1,2]);   offsetAvolt = BiasRctd.DcDiffOffsets.E12AVolt(iCalibTimeH);
          elseif isequal(antennas, [1,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E13AVolt(iCalibTimeH);
          elseif isequal(antennas, [2,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E23AVolt(iCalibTimeH);
          else
            error('BICAS:Assertion:IllegalArgument', 'Illegal argument "ssid".');
          end

        case 'AC_DIFF'

          if     isAchg == 0
            itfAvpiv    = BiasRctd.ItfSet.aclgAvpiv{iCalibTimeL};
            kFtfIvpav   = obj.BiasScalarGain.gammaIvpav.aclg;
            offsetAvolt = 0;
          elseif isAchg == 1
            itfAvpiv    = BiasRctd.ItfSet.achgAvpiv{iCalibTimeL};
            kFtfIvpav   = obj.BiasScalarGain.gammaIvpav.achg;
            offsetAvolt = 0;
          elseif isnan(isAchg)
            % CASE: AC GAIN unknown when it is NEEDED (i.e. when AC data).
            itfAvpiv    = obj.NAN_TF;
            kFtfIvpav   = NaN;
            offsetAvolt = NaN;
          else
            error('BICAS:Assertion:IllegalArgument', ...
              'Illegal argument isAchg=%g.', isAchg)
          end

        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            ['Illegal argument "ssid" implies illegal ASID category="%s".', ...
            ' Can not obtain calibration data for this type of signal.'], ...
            asidCategory)
      end

      %============================================
      % Modify return values for specific settings
      %============================================
      if obj.biasOffsetsDisabled && ~isnan(offsetAvolt)
        % NOTE: Overwrites "offsetAvolt".
        offsetAvolt = 0;
      end
      if obj.useBiasTfScalar
        % NOTE: Overwrites "itfAvpiv".
        itfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
      end
    end



    % Return LFR ITF, but also handle the case that should never happen for
    % actual non-NaN data (LSF F3 + BLTS 4 or 5) and return an ITF that only
    % returns NaN instead. BICAS may still iterate over that combination when
    % calibrating.
    %
    % RETURN VALUE
    % ============
    % itfIvpt
    %       LFR ITF, TM-->ivolt
    %
    function itfIvpt = get_LFR_ITF(obj, iBlts, iLsf, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      bicas.proc.L1L2.cal.utils.assert_iLsf( iLsf)
      assert(isscalar(iNonBiasRct))
      assert(iNonBiasRct >= 1, 'Illegal iNonBiasRct=%g', iNonBiasRct)
      % No assertion on zvcti2 unless used (determined later).

      %==================================================
      % The only place to potentially make use of zvcti2
      %==================================================
      if obj.useZvcti2
        % ASSERTIONS
        assert(isscalar(zvcti2), ...
          'BICAS:IllegalArgument:Assertion', ...
          'Argument zvcti2 is not scalar.')
        assert(zvcti2 >= 0, ...
          'BICAS:IllegalArgument:Assertion', ...
          ['Illegal argument zvcti2=%g', ...
          ' (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))'], ...
          zvcti2)
        assert(iLsf == zvcti2+1, ...
          'BICAS:IllegalArgument:Assertion', ...
          'zvcti2+1=%i != iLsf=%i (before overwriting iLsf)', ...
          zvcti2+1, iLsf)

        % NOTE: Override earlier iLsf.
        % NOTE: This is the only place zvcti2 is used in this class.
        iLsf = zvcti2 + 1;
      end



      if (iLsf == 4) && ismember(iBlts, [4,5])
        % CASE: F3 and BLTS={4,5}

        % IMPLEMENTATION NOTE: There is no tabulated LFR TF for this case and
        % the h/w does not support it, so an accurate TF can not be returned
        % even in principle. However, the BICAS implementation (2025-07-11)
        % iterates over all 5x BLTS channels no matter the value of LSF, e.g.
        % for F3 (iLsf=4) when there is only real data on BLTS1-3. It can
        % therefore request "calibration" for this case anyway, even if it means
        % calibrating an empty channel (converting NaN values to NaN). For this
        % reason, the code can not raise an exception for this case.
        itfIvpt = obj.NAN_TF;
      else
        LfrRctdCa = obj.Rctdc.get_RCTD_CA('LFR');

        % ASSERTION
        % IMPLEMENTATION NOTE: Anonymous function below will fail at a
        % later stage if these assertions are false. Checking for these
        % criteria here makes it easier to understand these particular
        % types of error.
        assert(iNonBiasRct <= numel(LfrRctdCa), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa is too small for argument iLfrRctd=%g.', ...
          ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
          ' value is larger than glob. attr. CALIBRATION TABLE allows.'], ...
          iNonBiasRct)
        assert(~isempty(LfrRctdCa{iNonBiasRct}), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa contains no RCT data corresponding', ...
          ' to argument iNonBiasRct=%g. This may indicate that', ...
          ' a zVar CALIBRATION_TABLE_INDEX(:,1) value is wrong or', ...
          ' that BICAS did not try to load the corresponding RCT', ...
          ' in glob. attr. CALIBRATION_TABLE.'], ...
          iNonBiasRct)

        itfIvpt = LfrRctdCa{iNonBiasRct}.ItfModifIvptCaCa{iLsf}{iBlts};
      end
    end    % get_LFR_ITF()



    % ARGUMENTS AND RETURN VALUES
    % ===========================
    % See get_LFR_ITF().
    %
    function itfIvpt = get_TDS_CWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(iNonBiasRct >= 1)
      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-CWF calibration never uses ZVCTI2.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        RctdCa        = obj.Rctdc.get_RCTD_CA('TDS-CWF');
        tdsFactorIvpt = RctdCa{iNonBiasRct}.factorsIvpt(iBlts);

        itfIvpt = @(omegaRps) (tdsFactorIvpt * ones(size(omegaRps)));
      else
        % CASE: BLTS 4-5 which TDS does NOT support (forbidden in h/w).

        itfIvpt = obj.NAN_TF;
      end
    end    % get_TDS_CWF_ITF()



    % ARGUMENTS AND RETURN VALUES
    % ===========================
    % See get_LFR_ITF().
    %
    function itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(iNonBiasRct >= 1)
      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-RSWF calibration never uses ZVCTI2.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        %======================================
        % Create combined ITF for TDS and BIAS
        %======================================
        RctdCa     = obj.Rctdc.get_RCTD_CA('TDS-RSWF');
        itfIvpt = RctdCa{iNonBiasRct}.itfModifIvptCa{iBlts};

      else
        % CASE: BLTS 4-5 which TDS does not support (forbidden in h/w).

        itfIvpt = obj.NAN_TF;
      end
    end    % get_TDS_RSWF_ITF()



    function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalibL = bicas.proc.L1L2.cal.utils.get_calibration_time_index(...
        Epoch, BiasRctdCa{1}.epochL);
    end



    function iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalibH = bicas.proc.L1L2.cal.utils.get_calibration_time_index(...
        Epoch, BiasRctdCa{1}.epochH);
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
