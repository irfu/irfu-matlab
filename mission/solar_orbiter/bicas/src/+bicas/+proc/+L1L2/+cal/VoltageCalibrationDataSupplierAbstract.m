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



    % Return BIAS ITF, BIAS scalar (ITF) calibration factor, and BIAS offset.
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
    % RETURN VALUE
    % ============
    % kItfAvpiv
    %       Multiplication factor "k" that represents/replaces the ITF for
    %       scalar calibration.
    % offsetAvolt
    %       Offset intended to be used with ITF, but could be used with scalar
    %       calibration too.
    %
    [itfAvpiv, kItfAvpiv, offsetAvolt] = get_BIAS_ITF_and_offset(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)

    itfIvpt = get_LFR_ITF(     obj, iBlts, NbriFpa, NbciFpa, iLsf)
    itfIvpt = get_TDS_CWF_ITF( obj, iBlts, NbriFpa, NbciFpa)
    itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, NbriFpa, NbciFpa)

    iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
    iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)



  end    % methods(Access=public)



end
