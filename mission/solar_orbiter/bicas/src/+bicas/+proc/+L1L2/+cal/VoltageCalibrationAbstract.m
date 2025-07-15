%
% Abstract superclass to bicas.proc.L1L2.cal.VoltageCalibrationImpl to (at
% least) facilitate mocking in automated tests.
%
% Class for functions which calibrate voltages.
%
% NOTE: Might not implement all abstract methods for all subclass methods which
% actually used externally. This needs to be done if automated tests (which
% should use another subclass) require access to more methods.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) VoltageCalibrationAbstract



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract, Access=public)



    bltsSamplesAVoltCa = calibrate_voltage_TM_to_avolt(obj, ...
        dtSec, bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, zvcti, ufv)



    % IMPLEMENTATION NOTE: This method technically has the same purpose as
    % bicas.proc.L1L2.cal.CurrentCalibrationAbstract.get_BIAS_calibration_time_index_L().
    iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
    iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)



  end    % methods(Access=public)



end
