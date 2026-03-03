%
% Subclass to be used by automated tests.
%
% Meant to be extended as more functionality is needed for automated tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VoltageCalibrationDataSupplierTest < bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierAbstract



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % BIAS
    itfBiasAvpiv
    kItfBiasAvpiv
    offsetAvolt

    % LFR/TDS
    itfLfrAvpiv
    itfTdsCwfAvpiv
    itfTdsRswfAvpiv
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = VoltageCalibrationDataSupplierTest(A)
      arguments
        % NOTE: All default vaulues are NaN to make it more likely that bad
        % tests fail.
        A.itfBiasAvpiv    = bicas.const.NAN_TF;
        A.kItfBiasAvpiv   = NaN;
        A.offsetAvolt     = NaN;
        A.itfLfrAvpiv     = bicas.const.NAN_TF;
        A.itfTdsCwfAvpiv  = bicas.const.NAN_TF;
        A.itfTdsRswfAvpiv = bicas.const.NAN_TF;
      end

      obj.itfBiasAvpiv    = A.itfBiasAvpiv;
      obj.kItfBiasAvpiv   = A.kItfBiasAvpiv;
      obj.offsetAvolt     = A.offsetAvolt;
      obj.itfLfrAvpiv     = A.itfLfrAvpiv;
      obj.itfTdsCwfAvpiv  = A.itfTdsCwfAvpiv;
      obj.itfTdsRswfAvpiv = A.itfTdsRswfAvpiv;
    end



    function [itfAvpiv, kItfAvpiv, offsetAvolt] = get_BIAS_ITF_and_offset(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)

      itfAvpiv    = obj.itfBiasAvpiv;
      kItfAvpiv   = obj.kItfBiasAvpiv;
      offsetAvolt = obj.offsetAvolt;
    end



    function itfIvpt = get_LFR_ITF(obj, iBlts, NbriFpa, NbciFpa, iLsf)
      assert(isa(NbriFpa, "bicas.utils.FPArray") & isscalar(NbriFpa))
      assert(isa(NbciFpa, "bicas.utils.FPArray") & isscalar(NbciFpa))

      itfIvpt = obj.itfLfrAvpiv;
    end



    function itfIvpt = get_TDS_CWF_ITF(obj, iBlts, NbriFpa, NbciFpa)
      assert(isa(NbriFpa, "bicas.utils.FPArray") & isscalar(NbriFpa))
      assert(isa(NbciFpa, "bicas.utils.FPArray") & isscalar(NbciFpa))

      itfIvpt = obj.itfTdsCwfAvpiv;
    end



    function itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, NbriFpa, NbciFpa)
      assert(isa(NbriFpa, "bicas.utils.FPArray") & isscalar(NbriFpa))
      assert(isa(NbciFpa, "bicas.utils.FPArray") & isscalar(NbciFpa))

      itfIvpt = obj.itfTdsRswfAvpiv;
    end



    function iCalibL = get_BIAS_calibration_time_index_L(obj, tt2000)
      iCalibL = ones(size(tt2000));
    end



    function iCalibH = get_BIAS_calibration_time_index_H(obj, tt2000)
      iCalibH = ones(size(tt2000));
    end



  end    % methods(Access=public)



end
