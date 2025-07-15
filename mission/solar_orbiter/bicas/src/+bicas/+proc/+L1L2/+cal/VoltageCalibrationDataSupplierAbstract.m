%
% Abstract superclass to bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierImpl
% to (at least) facilitate mocking in automated tests.
%
% Returns voltage calibration data without using it.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) VoltageCalibrationDataSupplierAbstract



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract, Access=public)



    [itfAvpiv, offsetAvolt] = get_BIAS_ITF_and_offset(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)

    itfIvpt = get_LFR_ITF(obj, iBlts, iLsf, iNonBiasRct, zvcti2)
    itfIvpt = get_TDS_CWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)
    itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)

    iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
    iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)



  end    % methods(Access=public)



end
