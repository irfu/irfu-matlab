%
% Class for functions that calibrate data.
% An instance of this class contains
%   (1) relevant settings (loaded from BSO) on how to calibrate data, and
%   (2) calibration data.
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
% Does not implement parasitic capacitance due to lack of calibration values (at
% least). Should not need to implement according to Thomas Chust(?).
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
% (4) to detect changing RCT formats, in particular in RCTs from non-BIAS teams.
% --
% NOTE: All calibration functions of measured data are assumed to accept data
% from all BLTS (1-5), i.e. including TDS, in order to reduce the number
% assumptions that the calling code needs to make.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% Note: See readme.txt.
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
% CALIBRATION_TABLE{CALIBRATION_TABLE_INDEX{iRecord, 1} + 1}
%     == RCT filename
% CALIBRATION_TABLE_INDEX{iRecord, 2}
%     == ZVCTI2
%     == Index/pointer to some calibration value(s) to use in the corresponding
%        RCT. The exact interpretation depends on the RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef Cal < handle
  % All methods 2024-09-05
  % ----------------------
  % function obj = Cal(...
  % function biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
  % function biasCurrentAAmpere = calibrate_current_HK_TM_to_aampere(obj, ...
  % function bltsSamplesAVoltCa = calibrate_voltage_all(obj, ...
  %   Delegates to calibrate_voltage_*() but also handles
  %   allVoltageCalibDisabled, ufv, useGact(=false) itself.
  % function bltsSamplesAVoltCa = calibrate_voltage_BIAS_LFR(obj, ...
  % function bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_CWF(obj, ...
  % function bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_RSWF(obj, ...
  % function iCalib = get_BIAS_calibration_time_L(obj, Epoch)
  % function iCalib = get_BIAS_calibration_time_H(obj, Epoch)
  % function BiasCalibData = get_BIAS_calib_data(obj, ...
  % function lfrItfIvpt = get_LFR_ITF(obj, iLfrRctd, iBlts, iLsf)
  % function [CalData] = get_BIAS_LFR_calib_data(obj, CalSettings, iNonBiasRct, zvcti2)
  % function biasCurrentTm = calibrate_current_sampere_to_TM(currentSAmpere)
  %
  %
  %
  % PROPOSAL: Separate classes for different types of data. At least separate
  %           for LFR, TDS-CWF, TDS-RSWF.
  %   NOTE: Should then also have separate class for current calibration.
  %   PROPOSAL: Functions
  %     bltsSamplesAVoltCa = calibrate_voltage_BIAS_LFR(obj, ...
  %     bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_CWF(obj, ...
  %     bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_RSWF(obj, ...
  %   should then be separate implementations of abstract superclass mathod.
  %     PROPOSAL: Separate non-abstract superclass method (which handles
  %       allVoltageCalibDisabled, ufv, useGact(=false)) calls the subclass
  %        method.
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
  %                   bicas.proc.L1L2.cal.Cal encloses/contains/"secretly uses"
  %                   bicas.tf.apply_TF_freq().
  %             (2) For testing, one wants to verify the path (both ways) between
  %                   bicas.proc.L1L2.cal.Cal and
  %                   bicas.tf.apply_TF_freq().
  %             ==> One wants to test one unit of code at a time, but what a
  %                "unit" is ambiguous:
  %                 one wants small units of code
  %                 unit is ambiguous when a unit uses/call other unit(s).
  %
  %
  %
  % PROPOSAL: Move (charge) current calibration to separate class.
  %   NOTE: Functions
  %       calibrate_current_TM_to_aampere()
  %           Uses BiasRctd == Uses RCT.
  %       calibrate_current_HK_TM_to_aampere()
  %           Uses
  %               obj.HkBiasCurrent.gainAapt
  %               obj.HkBiasCurrent.offsetTm
  %       calibrate_current_sampere_to_TM()
  %           Static
  %           Uses solo.hwzv.const.TM_PER_SAMPERE
  %   PRO: Class is large, ~1200 rows.
  %   PRO: Remaining class becomes entirely about voltage calibration.
  %   CON: Use needs to instantiate two calibration objects.
  %   TODO-DEC: Name of new class?
  %       ~cal_curr
  %   PROPOSAL: Rename remaining class: Only about voltage calibration.
  %       ~cal_volt
  %
  %
  %
  % TODO-DEC: How distribute the calibration formulas/algorithms between
  %   (1) calibrate_* functions,
  %   (2) functions that select relevant calibration data (get_BIAS_calib_data)
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
  % TODO-NI: Where does the parasitic capacitance TF fit into the calibration formulas?
  % TODO-NI: What parasitic capacitance value(s) should one use?
  % PROPOSAL: Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
  %
  % PROPOSAL: Store all versions of TFs internally.
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
  % PROPOSAL: Refactor to use a struct constant for those arguments to
  %           bicas.tf.apply_TF() which are constant.
  %
  % BUG: Can likely not handle data with SSID = Unknown or 2.5V Ref, at least
  %      not for LFR.
  %   PROPOSAL: Tests.



  properties(Access=private, Constant)

    % Local TF constant for convenience.
    NAN_TF = @(omegaRps) (omegaRps * NaN);

  end



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
    HkBiasCurrent



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



  %###########################################################################



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
    %    submitting it to bicas.proc.L1L2.cal.Cal
    % ** it simplifies the constructor.
    %
    function obj = Cal(Rctdc, useGactRct, useZvcti2, Bso)

      % ASSERTIONS: Arguments
      assert(isscalar(useZvcti2))
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
      %       (more readable code).
      %   ** Potentially gives faster access to values (better
      %      performance).
      %==================================================================
      obj.HkBiasCurrent.offsetTm         = Bso.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM');
      obj.HkBiasCurrent.gainAapt         = Bso.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_AAPT');

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

      obj.allVoltageCalibDisabled        = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLE');
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



    % Convert/calibrate TC bias current: TM units --> physical units.
    %
    % NOTE: This is the normal way of obtaining bias current in physical
    % units (as opposed to HK bias current).
    %
    % ARGUMENTS
    % =========
    % iCalibTimeL
    %       Has to be same size as "biasCurrentTm".
    %
    function biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna, iCalibTimeL)

      assert(isscalar(iAntenna))
      assert(isequaln(...
        size(biasCurrentTm), ...
        size(iCalibTimeL)))

      %==============================
      % Obtain calibration constants
      %==============================
      BiasRctdCa    = obj.Rctdc.get_RCTD_CA('BIAS');
      offsetAAmpere = BiasRctdCa{1}.Current.offsetsAAmpere(iCalibTimeL, iAntenna);
      gainAapt      = BiasRctdCa{1}.Current.gainsAapt(     iCalibTimeL, iAntenna);

      % CALIBRATE
      %
      % LINEAR FUNCTION
      biasCurrentAAmpere = offsetAAmpere + gainAapt .* double(biasCurrentTm);
    end



    % Convert/calibrate diagnostic HK TM bias current values to physical
    % units. Refers to BIAS HK ZVs HK_BIA_BIAS1/2/3.
    %
    %
    % NOTES
    % =====
    % IMPORTANT NOTE: The HK bias current values are measured onboard but
    % are only meant as DIAGNOSTIC values, NOT AS THE PROPER BIAS CURRENT
    % values for nominal use. Therefore the values should only be seen as
    % approximate.
    %
    % NOTE: Walter Puccio, IRF-U 2019-09-06: Values are measured on the
    % order of once per second (and sent back as HK even more rarely).
    % Expect errors on the order of 5%.
    %
    % NOTE: The calibration data are NOT stored in the BIAS RCT.
    %
    % NOTE: The conversion function can be found in the BIAS specification,
    % sections 3.4.4.{1-3} ("BIAS1" etc) under "Telemetry". (Not to be
    % confused with the corresponding telecommands.). The conversion
    % functions are identical for all three probes.
    %
    function biasCurrentAAmpere = calibrate_current_HK_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna)

      % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK.
      % Not strictly required, but the variable has to be some integer.
      assert(isa(biasCurrentTm, 'uint16'))

      %=============================================================
      % CALIBRATE
      % ---------
      % Unsigned integer which represents ~signed integer.
      % ==> Intervals 0..0x7FFF and 0x8000...0xFFFF need to
      %     "change places".
      % ==> Need to flip bit representing sign to have one interval
      %     0...0xFFFF with monotonic function for TM-to-calibrated
      %     values.
      %=============================================================
      biasCurrentTm      = bitxor(biasCurrentTm, hex2dec('8000'));    % FLIP BIT
      biasCurrentAAmpere = obj.HkBiasCurrent.gainAapt(iAntenna) * ...
        (biasCurrentTm + obj.HkBiasCurrent.offsetTm(iAntenna));    % LINEAR FUNCTION
    end



    % Calibrate all voltages. Function will choose the more specific
    % algorithm internally.
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
        dtSec, bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, ...
        zvcti, ufv)

      % ASSERTIONS
      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))
      assert(iscell(bltsSamplesTmCa))
      assert(isvector(bltsSamplesTmCa))
      irf.assert.sizes(zvcti, [1,2])
      assert(islogical(ufv) && isscalar(ufv))



      % Set iNonBiasRct, zvcti2 by extracting values from zvcti or emulating it.
      if obj.useGactRct
        % NOTE: Incrementing by one (index into MATLAB array).
        iNonBiasRct = 1 + zvcti(1,1);
      else
        iNonBiasRct = 1;
      end
      % NOTE: NOT incrementing value by one, since the variable's meaning
      % can vary between LFR, TDS-CWF, TDS-RSWF.
      zvcti2 = zvcti(1,2);



      if obj.allVoltageCalibDisabled || ufv

        bltsSamplesAVoltCa = cell(size(bltsSamplesTmCa));

        for i = 1:numel(bltsSamplesTmCa)
          if obj.allVoltageCalibDisabled
            % CASE: Set voltages to TM values.
            bltsSamplesAVoltCa{i} = double(bltsSamplesTmCa{i});
          end
          if ufv
            % CASE: Set voltages to NaN.

            % IMPLEMENTATION NOTE: Potentially overwrites TM value set in above
            % "if" statement.
            bltsSamplesAVoltCa{i} = nan(size(bltsSamplesTmCa{i}));
          end
        end

      else

        if isLfr
          %===========
          % CASE: LFR
          %===========
          bltsSamplesAVoltCa = obj.calibrate_voltage_BIAS_LFR(...
            dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2);
        else
          %===========
          % CASE: TDS
          %===========
          if isTdsCwf
            % CASE: TDS CWF
            bltsSamplesAVoltCa = obj.calibrate_voltage_BIAS_TDS_CWF(...
              dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2);
          else
            % CASE: TDS RSWF
            bltsSamplesAVoltCa = obj.calibrate_voltage_BIAS_TDS_RSWF(...
              dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2);
          end
        end

      end
    end



    % ARGUMENTS
    % =========
    % samplesTm    : 1D cell array of numeric 1D arrays.
    % samplesAVolt : 1D cell array of numeric 1D arrays.
    %
    function bltsSamplesAVoltCa = calibrate_voltage_BIAS_LFR(obj, ...
        dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2)

      % ASSERTIONS
      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))
      % IMPLEMENTATION NOTE: bicas.proc.L1L2.CalibrationSettings permits TDS
      % data for which iLsf=NaN.
      bicas.proc.L1L2.cal.utils.assert_iLsf(CalSettings.iLsf)
      irf.assert.vector(bltsSamplesTmCa)
      assert(iscell(bltsSamplesTmCa) && isvector(bltsSamplesTmCa))
      irf.assert.vector(dtSec)
      assert(numel(bltsSamplesTmCa) == numel(dtSec))



      %=============================
      % Obtain all calibration data
      %=============================
      CalibData = obj.get_BIAS_LFR_calib_data(...
        CalSettings, iNonBiasRct, zvcti2);

      %====================================
      % CALIBRATE: LFR TM --> TM --> avolt
      %====================================
      bltsSamplesAVoltCa = cell(size(bltsSamplesTmCa));
      for i = 1:numel(bltsSamplesTmCa)

        % APPLY TRANSFER FUNCTION (BIAS + LFR)
        tempSamplesAVolt = bicas.tf.apply_TF(...
          dtSec(i), ...
          bltsSamplesTmCa{i}(:), ...
          CalibData.itfAvpt, ...
          'method',                  obj.tfMethod, ...
          'detrendingDegreeOf',      CalibData.detrendingDegreeOf, ...
          'retrendingEnabled',       CalibData.retrendingEnabled, ...
          'tfHighFreqLimitFraction', CalibData.itfHighFreqLimitFraction, ...
          'kernelEdgePolicy',        obj.kernelEdgePolicy, ...
          'kernelHannWindow',        obj.kernelHannWindow, ...
          'snfEnabled',              obj.snfEnabled, ...
          'snfSubseqMinSamples',     obj.snfSubseqMinSamples);

        % ADD BIAS offset
        bltsSamplesAVoltCa{i} = tempSamplesAVolt + CalibData.BiasCalibData.offsetAVolt;
      end
    end



    % ARGUMENTS
    % =========
    % See calibrate_voltage_BIAS_LFR.
    %
    function bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_CWF(obj, ...
        dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2)

      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))

      iBlts        = CalSettings.iBlts;
      Ssid         = CalSettings.Ssid;
      isAchg       = CalSettings.isAchg;
      iCalibTimeL  = CalSettings.iCalibTimeL;
      iCalibTimeH  = CalSettings.iCalibTimeH;

      % ASSERTIONS
      irf.assert.vector(dtSec)
      assert(iscell(bltsSamplesTmCa) && isvector(bltsSamplesTmCa))
      assert(numel(bltsSamplesTmCa) == numel(dtSec))
      assert(iNonBiasRct >= 1)

      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-CWF calibration never uses ZVCTI2.')
      end

      % Initialize empty output variable.
      bltsSamplesAVoltCa = cell(size(bltsSamplesTmCa));

      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        %==============================
        % Obtain calibration constants
        %==============================
        % NOTE: AC low/high gain is irrelevant for TDS. Argument value is
        % therefore arbitrary.
        BiasCalibData = obj.get_BIAS_calib_data(...
          Ssid, isAchg, iCalibTimeL, iCalibTimeH);

        if obj.lfrTdsTfDisabled
          tdsFactorIvpt = 1;
        else
          RctdCa        = obj.Rctdc.get_RCTD_CA('TDS-CWF');
          tdsFactorIvpt = RctdCa{iNonBiasRct}.factorsIvpt(iBlts);
        end

        for i = 1:numel(bltsSamplesTmCa)

          %===============================================
          % CALIBRATE: TDS TM --> TDS/BIAS interface volt
          %===============================================
          % MULTIPLICATION
          tempSamplesIVolt = tdsFactorIvpt * bltsSamplesTmCa{i};

          %=====================================================
          % CALIBRATE: TDS/BIAS interface volt --> antenna volt
          %=====================================================
          % APPLY TRANSFER FUNCTION (for BIAS, but not for TDS-CWF)
          tempSamplesAVolt = bicas.tf.apply_TF(...
            dtSec(i), ...
            tempSamplesIVolt(:), ...
            BiasCalibData.itfAvpiv, ...
            'method',                  obj.tfMethod, ...
            'detrendingDegreeOf',      obj.dcDetrendingDegreeOf, ...
            'retrendingEnabled',       obj.dcRetrendingEnabled, ...
            'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction, ...
            'kernelEdgePolicy',        obj.kernelEdgePolicy, ...
            'kernelHannWindow',        obj.kernelHannWindow, ...
            'snfEnabled',              obj.snfEnabled, ...
            'snfSubseqMinSamples',     obj.snfSubseqMinSamples);

          % ADD BIAS OFFSET
          bltsSamplesAVoltCa{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
        end

      else
        % CASE: BLTS 4-5 which TDS does NOT support.

        for i = 1:numel(bltsSamplesTmCa)
          % Always return NaN.
          bltsSamplesAVoltCa{i} = NaN * bltsSamplesTmCa{i};
        end
      end

    end



    % ARGUMENTS
    % =========
    % See calibrate_voltage_BIAS_LFR.
    %
    function bltsSamplesAVoltCa = calibrate_voltage_BIAS_TDS_RSWF(obj, ...
        dtSec, bltsSamplesTmCa, CalSettings, iNonBiasRct, zvcti2)

      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))

      iBlts        = CalSettings.iBlts;
      Ssid         = CalSettings.Ssid;
      isAchg       = CalSettings.isAchg;
      iCalibTimeL  = CalSettings.iCalibTimeL;
      iCalibTimeH  = CalSettings.iCalibTimeH;

      % ASSERTIONS
      irf.assert.vector(dtSec)
      assert(iscell(bltsSamplesTmCa) && isvector(bltsSamplesTmCa))
      assert(numel(bltsSamplesTmCa) == numel(dtSec))
      assert(iNonBiasRct >= 1)

      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-RSWF calibration never uses ZVCTI2.')
      end

      %==============================
      % Obtain calibration constants
      %==============================
      % NOTE: Low/high gain is irrelevant for TDS. Argument value
      % arbitrary.
      BiasCalibData = obj.get_BIAS_calib_data(...
        Ssid, isAchg, iCalibTimeL, iCalibTimeH);

      % Initialize empty output variable.
      bltsSamplesAVoltCa = cell(size(bltsSamplesTmCa));
      if ismember(iBlts, [1,2,3])

        %======================================
        % Create combined ITF for TDS and BIAS
        %======================================
        if obj.lfrTdsTfDisabled
          tdsItfIvpt = @(omegaRps) (ones(omegaRps));
        else
          RctdCa     = obj.Rctdc.get_RCTD_CA('TDS-RSWF');
          tdsItfIvpt = RctdCa{iNonBiasRct}.itfModifIvptCa{iBlts};
        end

        itfAvpt = @(omegaRps) (...
          tdsItfIvpt(omegaRps) ...
          .* ...
          BiasCalibData.itfAvpiv(omegaRps));

        %====================================
        % CALIBRATE: TDS TM --> antenna volt
        %====================================
        % APPLY TRANSFER FUNCTION (BIAS + TDS-RSWF)
        for i = 1:numel(bltsSamplesTmCa)
          tempSamplesAVolt = bicas.tf.apply_TF(...
            dtSec(i), ...
            bltsSamplesTmCa{i}(:), ...
            itfAvpt, ...
            'method',                  obj.tfMethod, ...
            'detrendingDegreeOf',      obj.dcDetrendingDegreeOf, ...
            'retrendingEnabled',       obj.dcRetrendingEnabled, ...
            'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction, ...
            'kernelEdgePolicy',        obj.kernelEdgePolicy, ...
            'kernelHannWindow',        obj.kernelHannWindow, ...
            'snfEnabled',              obj.snfEnabled, ...
            'snfSubseqMinSamples',     obj.snfSubseqMinSamples);

          % ADD BIAS OFFSET
          bltsSamplesAVoltCa{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
        end
      else
        for i = 1:numel(bltsSamplesTmCa)
          % CASE: BLTS 4-5 which TDS does not support.
          % Always return NaN.
          bltsSamplesAVoltCa{i} = NaN * bltsSamplesTmCa{i};
        end
      end

    end



    function iCalib = get_BIAS_calibration_time_L(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalib = bicas.proc.L1L2.cal.utils.get_calibration_time(...
        Epoch, BiasRctdCa{1}.epochL);
    end



    function iCalib = get_BIAS_calibration_time_H(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalib = bicas.proc.L1L2.cal.utils.get_calibration_time(...
        Epoch, BiasRctdCa{1}.epochH);
    end



    % Return subset of already loaded BIAS calibration data, for specified
    % settings.
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
    function BiasCalibData = get_BIAS_calib_data(obj, ...
        Ssid, isAchg, iCalibTimeL, iCalibTimeH)

      % PROPOSAL: Log warning message when simultaneously isAchg=NaN
      % and the value is needed.

      % ASSERTION
      assert(bicas.sconst.is_SSID(Ssid) & isscalar(Ssid))
      assert(bicas.sconst.SSID_is_ASR(Ssid))
      assert(isscalar(isAchg) && isnumeric(isAchg))
      assert(isscalar(iCalibTimeL))
      assert(isscalar(iCalibTimeH))

      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');
      BiasRctd   = BiasRctdCa{1};

      %###################################################################
      % kIvpav = Multiplication factor "k" that represents/replaces the
      % (forward) transfer function.
      asid         = bicas.sconst.SSID_ASR_to_ASID(Ssid);
      asidCategory = bicas.sconst.get_ASID_category(asid);
      antennas     = bicas.sconst.get_ASID_antennas(asid);
      switch(asidCategory)
        case 'DC_SINGLE'

          % NOTE: List of ITFs for different times.
          biasItfAvpiv = BiasRctd.ItfSet.dcSingleAvpiv{iCalibTimeL};
          kFtfIvpav    = obj.BiasScalarGain.alphaIvpav;
          offsetAVolt  = BiasRctd.dcSingleOffsetsAVolt(iCalibTimeH, antennas);

        case 'DC_DIFF'

          biasItfAvpiv = BiasRctd.ItfSet.dcDiffAvpiv{iCalibTimeL};
          kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
          if     isequal(antennas, [1,2]);   offsetAVolt = BiasRctd.DcDiffOffsets.E12AVolt(iCalibTimeH);
          elseif isequal(antennas, [1,3]);   offsetAVolt = BiasRctd.DcDiffOffsets.E13AVolt(iCalibTimeH);
          elseif isequal(antennas, [2,3]);   offsetAVolt = BiasRctd.DcDiffOffsets.E23AVolt(iCalibTimeH);
          else
            error('BICAS:Assertion:IllegalArgument', ...
              'Illegal Ssid.');
          end

        case 'AC_DIFF'

          if     isAchg == 0
            biasItfAvpiv = BiasRctd.ItfSet.aclgAvpiv{iCalibTimeL};
            kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.aclg;
            offsetAVolt  = 0;
          elseif isAchg == 1
            biasItfAvpiv = BiasRctd.ItfSet.achgAvpiv{iCalibTimeL};
            kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.achg;
            offsetAVolt  = 0;
          elseif isnan(isAchg)
            % CASE: GAIN unknown when it is NEEDED for calibration.
            biasItfAvpiv = bicas.proc.L1L2.cal.Cal.NAN_TF;
            kFtfIvpav    = NaN;
            offsetAVolt  = NaN;
          else
            error('BICAS:Assertion:IllegalArgument', ...
              'Illegal argument isAchg=%g.', isAchg)
          end

        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            ['Illegal argument Ssid.Asid.category=%s.', ...
            ' Can not obtain calibration data for this type of signal.'], ...
            asidCategory)
      end

      if obj.biasOffsetsDisabled && ~isnan(offsetAVolt)
        % NOTE: Overwrites "offsetAVolt".
        offsetAVolt = 0;
      end
      if obj.useBiasTfScalar
        % NOTE: Overwrites "biasItfAvpiv".
        biasItfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
      end
      %###################################################################



      BiasCalibData.itfAvpiv    = biasItfAvpiv;
      BiasCalibData.offsetAVolt = offsetAVolt;

    end



    % Obtain LFR ITF, but handle the case that should never happen for
    % actual non-NaN data (LSF F3 + BLTS 4 or 5) and return a TF that only
    % returns NaN instead. BICAS may still iterate over that combination
    % though when calibrating.
    %
    function lfrItfIvpt = get_LFR_ITF(obj, iLfrRctd, iBlts, iLsf)
      % ASSERTIONS
      assert(iLfrRctd >= 1)
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      bicas.proc.L1L2.cal.utils.assert_iLsf(iLsf)

      if (iLsf == 4) && ismember(iBlts, [4,5])
        % CASE: F3 and BLTS={4,5}

        % NOTE: There is no tabulated LFR TF and no such combination
        % signal route, so the TF can not be returned even in principle.
        lfrItfIvpt = bicas.proc.L1L2.cal.Cal.NAN_TF;
      else
        LfrRctdCa = obj.Rctdc.get_RCTD_CA('LFR');

        % ASSERTION
        % IMPLEMENTATION NOTE: Anonymous function below will fail at a
        % later stage if these assertions are false. Checking for these
        % criteria here makes it easier to understand these particular
        % types of error.
        assert(numel(LfrRctdCa) <= iLfrRctd, ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa is too small for argument iLfrRctd=%g.', ...
          ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
          ' value is larger than glob. attr. CALIBRATION TABLE allows.'], ...
          iLfrRctd)
        assert(~isempty(LfrRctdCa{iLfrRctd}), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa contains no RCT data corresponding', ...
          ' to argument iLfrRctd=%g. This may indicate that', ...
          ' a zVar CALIBRATION_TABLE_INDEX(:,1) value is wrong or', ...
          ' that BICAS did not try to load the corresponding RCT', ...
          ' in glob. attr. CALIBRATION_TABLE.'], ...
          iLfrRctd)

        lfrItfIvpt = LfrRctdCa{iLfrRctd}.ItfModifIvptCaCa{iLsf}{iBlts};
      end
    end



    % Return calibration data for LFR+BIAS calibration.
    %
    % RATIONALE
    % =========
    % Method exists to
    % (1) simplify & clarify calibrate_voltage_BIAS_LFR(),
    % (2) be useful for non-BICAS code to inspect the calibration data used
    %     for a particular calibration case, in particular the combined
    %     (LFR+BIAS) transfer functions.
    %     NOTE: This is also the reason why this method is public.
    %
    % IMPLEMENTATION NOTE: Return one struct instead of multiple return
    % values to make sure that the caller does not confuse the return values
    % with each other.
    function [CalData] = get_BIAS_LFR_calib_data(obj, CalSettings, iNonBiasRct, zvcti2)

      assert(isa(CalSettings, 'bicas.proc.L1L2.CalibrationSettings'))

      iBlts        = CalSettings.iBlts;
      Ssid         = CalSettings.Ssid;
      isAchg       = CalSettings.isAchg;
      iCalibTimeL  = CalSettings.iCalibTimeL;
      iCalibTimeH  = CalSettings.iCalibTimeH;
      iLsf         = CalSettings.iLsf;

      % ASSERTIONS
      assert(bicas.sconst.SSID_is_ASR(Ssid))
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



      CalData = struct();

      %====================================================
      % Obtain settings for bicas.tf.apply_TF()
      %====================================================
      if bicas.sconst.SSID_is_AC(Ssid)
        % IMPLEMENTATION NOTE: DC is (optionally) detrended via
        % bicas.tf.apply_TF() in the sense of a linear fit
        % being removed, TF applied, and then added back. That same
        % algorithm, or at least adding back the fit, is by its nature
        % inappropriate for non-lowpass filters, i.e. for AC. (The fit
        % can not be scaled with the 0 Hz signal amplitude)
        CalData.detrendingDegreeOf = obj.acDetrendingDegreeOf;
        CalData.retrendingEnabled  = false;   % NOTE: HARDCODED SETTING.
      else
        CalData.detrendingDegreeOf = obj.dcDetrendingDegreeOf;
        CalData.retrendingEnabled  = obj.dcRetrendingEnabled;
      end
      CalData.itfHighFreqLimitFraction = obj.itfHighFreqLimitFraction;

      %==============================
      % Obtain BIAS calibration data
      %==============================
      CalData.BiasCalibData = obj.get_BIAS_calib_data(...
        Ssid, isAchg, iCalibTimeL, iCalibTimeH);

      %========================================
      % Obtain (official) LFR calibration data
      %========================================
      if obj.lfrTdsTfDisabled
        CalData.lfrItfIvpt = @(omegaRps) (ones(size(omegaRps)));
      else
        CalData.lfrItfIvpt = obj.get_LFR_ITF(iNonBiasRct, iBlts, iLsf);
      end

      %======================================
      % Create combined ITF for LFR and BIAS
      %======================================

      CalData.itfAvpt = bicas.proc.L1L2.cal.utils.create_LFR_BIAS_ITF(...
        CalData.lfrItfIvpt, ...
        CalData.BiasCalibData.itfAvpiv, ...
        bicas.sconst.SSID_is_AC(Ssid), ...
        obj.itfAcConstGainLowFreqRps);
    end



  end    % methods(Access=private)



  %###########################################################################



  methods(Static, Access=public)



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



    % Convert "set current" to TC/TM units.
    %
    function biasCurrentTm = calibrate_current_sampere_to_TM(currentSAmpere)

      % ASSERTION: Input values are within range.
      % NOTE: max(...) ignores NaN, unless that is the only value, which
      % then becomes the max value.
      [maxAbsSAmpere, iMax] = max(abs(currentSAmpere(:)));
      if ~(isnan(maxAbsSAmpere) || (maxAbsSAmpere <= solo.hwzv.const.MAX_ABS_SAMPERE))
        error('BICAS:Assertion:IllegalArgument', ...
          ['Argument currentSAmpere (unit: set current/ampere)', ...
          ' contains illegally large value(s).', ...
          ' Largest found value is %g.'], ...
          currentSAmpere(iMax))
      end

      biasCurrentTm = currentSAmpere * solo.hwzv.const.TM_PER_SAMPERE;
    end



  end    % methods(Static, Access=public)



end
