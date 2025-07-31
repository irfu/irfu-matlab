%
% Nominal implementation of superclass for calibration currents.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef CurrentCalibrationImpl < bicas.proc.L1L2.cal.CurrentCalibrationAbstract
  % PROPOSAL: Automatic test code.
  % PROPOSAL: Move HK constants BSO-->bicas.const.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    biasRctdEpochL
    Current
    HkBiasCurrent
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = CurrentCalibrationImpl(BiasRctd, Bso)
      % PROPOSAL: Replace argument BiasRctd --> offsetsAampere, gainsAapt,
      %           epochL
      %   PRO: Good for testing.
      %   CON: Harder to call constructor.
      assert(isa(BiasRctd, 'bicas.proc.L1L2.cal.rct.RctDataBias'))

      obj.biasRctdEpochL         = BiasRctd.epochL;
      obj.Current                = BiasRctd.Current;

      obj.HkBiasCurrent.offsetTm = Bso.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM');
      obj.HkBiasCurrent.gainAapt = Bso.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_AAPT');
    end



    % Convert "set current" to TC/TM units.
    %
    function biasCurrentTm = calibrate_current_sampere_to_TM(obj, currentSampere)

      % ASSERTION: Input values are within range.
      % NOTE: max(...) ignores NaN, unless that is the only value, which
      % then becomes the max value.
      [maxAbsSampere, iMax] = max(abs(currentSampere(:)));
      if ~(isnan(maxAbsSampere) || (maxAbsSampere <= solo.hwzv.const.MAX_ABS_SAMPERE))
        error('BICAS:Assertion:IllegalArgument', ...
          ['Argument currentSampere (unit: set current/ampere)', ...
          ' contains illegally large value(s).', ...
          ' Largest found value is %g.'], ...
          currentSampere(iMax))
      end

      biasCurrentTm = currentSampere * solo.hwzv.const.TM_PER_SAMPERE;
    end



    % Convert/calibrate TC bias current: TM units --> physical units.
    %
    % NOTE: This is the normal way of obtaining bias current in physical
    %       units (as opposed to HK bias current).
    %
    % ARGUMENTS
    % =========
    % iCalibTimeL
    %       Has to be same size as "biasCurrentTm".
    %
    function biasCurrentAampere = calibrate_current_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna, iCalibTimeL)

      assert(isscalar(iAntenna))
      assert(isequaln(...
        size(biasCurrentTm), ...
        size(iCalibTimeL)))

      %==============================
      % Obtain calibration constants
      %==============================
      offsetAampere = obj.Current.offsetsAampere(iCalibTimeL, iAntenna);
      gainAapt      = obj.Current.gainsAapt(     iCalibTimeL, iAntenna);

      % CALIBRATE
      %
      % LINEAR FUNCTION
      biasCurrentAampere = offsetAampere + gainAapt .* double(biasCurrentTm);
    end



    % Convert/calibrate diagnostic HK TM bias current values to physical
    % units. Refers to BIAS HK ZVs HK_BIA_BIAS1/2/3.
    %
    %
    % NOTES
    % =====
    % NOTE: This function is not used by BICAS. It exists because it could be
    % useful for manual work and is related.
    %
    % IMPORTANT NOTE: The HK bias current values are measured onboard but are
    % only meant as DIAGNOSTIC values, NOT AS THE PROPER BIAS CURRENT values
    % for nominal use. Therefore the values should only be seen as approximate.
    %
    % NOTE: Walter Puccio, IRF-U 2019-09-06: Values are measured on the order
    % of once per second (and sent back as HK even more rarely). Expect errors
    % on the order of 5%.
    %
    % NOTE: The calibration data are NOT stored in the BIAS RCT.
    %
    % NOTE: The conversion function can be found in the BIAS specification,
    % sections 3.4.4.{1-3} ("BIAS1" etc) under "Telemetry". (Not to be confused
    % with the corresponding telecommands.). The conversion functions are
    % identical for all three probes.
    %
    function biasCurrentAampere = calibrate_current_HK_TM_to_aampere(obj, ...
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
      biasCurrentAampere = obj.HkBiasCurrent.gainAapt(iAntenna) * ...
        (biasCurrentTm + obj.HkBiasCurrent.offsetTm(iAntenna));    % LINEAR FUNCTION
    end



    function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
      iCalibL = bicas.proc.L1L2.cal.utils.get_calibration_time_index(...
        Epoch, obj.biasRctdEpochL);
    end



  end    % methods(Access=public)



end
