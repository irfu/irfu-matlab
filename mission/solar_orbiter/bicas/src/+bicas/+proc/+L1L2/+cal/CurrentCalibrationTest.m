%
% Alternative implementation for the purpose of automated tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef CurrentCalibrationTest < bicas.proc.L1L2.cal.CurrentCalibrationAbstract



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    tmPerSampere
    iCalibLDict
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = CurrentCalibrationTest(tmPerSampere, iCalibLDict)
      obj.tmPerSampere = tmPerSampere;
      obj.iCalibLDict  = iCalibLDict;
    end



    function biasCurrentTm = calibrate_current_sampere_to_TM(obj, currentSampere)
      % PROPOSAL: Change?
      biasCurrentTm = obj.tmPerSampere * currentSampere;
    end



    function biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna, iCalibTimeL)

      biasCurrentAAmpere = biasCurrentTm .* iCalibTimeL + iAntenna;
    end



    function biasCurrentAAmpere = calibrate_current_HK_TM_to_aampere(obj, ...
        biasCurrentTm, iAntenna)

      error("NOT IMPLEMENTED")
    end



    % IMPLEMENTATION NOTE: This method technically has the same purpose as
    % bicas.proc.L1L2.cal.VoltageCalibration.get_BIAS_calibration_time_index_L().
    function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
      iCalibL = obj.iCalibLDict(Epoch);
    end



  end    % methods(Access=public)



end
