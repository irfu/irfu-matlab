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



    biasCurrentTm = calibrate_current_sampere_to_TM(obj, currentSampere)

    biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna, iCalibTimeL)

    biasCurrentAAmpere = calibrate_current_HK_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna)

    % IMPLEMENTATION NOTE: This method technically has the same purpose as
    % bicas.proc.L1L2.cal.VoltageCalibration.get_BIAS_calibration_time_index_L().
    iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)



  end    % methods(Abstract, Access=public)



end
