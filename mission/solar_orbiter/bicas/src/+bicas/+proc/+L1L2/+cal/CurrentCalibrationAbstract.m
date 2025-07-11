%
% Abstract superclass to bicas.proc.L1L2.cal.CurrentCalibrationImpl to (at
% least) facilitate mocking in automated tests.
%
% Class for functions which calibrate currents.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) CurrentCalibrationAbstract



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract, Access=public)



    biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna, iCalibTimeL)

    biasCurrentAAmpere = calibrate_current_HK_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna)

    % IMPLEMENTATION NOTE: This method technically has the same purpose as
    % bicas.proc.L1L2.cal.VoltageCalibrationAbstract.get_BIAS_calibration_time_index_L().
    iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)



  end    % methods(Abstract, Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



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



  end    % methods(Static)



end
