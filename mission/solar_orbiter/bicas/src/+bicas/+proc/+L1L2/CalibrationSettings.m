%
% Dumb container for some calibration settings. Effectively a replacement for a
% recurring set of function arguments.
%
% NOTE: Not all fields have valid values for all types of data.
%       Ex: isAchg, iLsf.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef CalibrationSettings



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=public)
      iBlts
      ssid
      isAchg
      iCalibTimeL
      iCalibTimeH
      iLsf
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = CalibrationSettings(...
        iBlts, ssid, isAchg, iCalibTimeL, iCalibTimeH, iLsf)

      % PROPOSAL: Assertions.
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(bicas.sconst.is_SSID(ssid) & isscalar(ssid))
      assert(isnan(isAchg) || ismember(isAchg, [0, 1]))
      if ~isnan(iLsf)
        % CASE: LFR data (not TDS)
        bicas.proc.L1L2.cal.utils.assert_iLsf(iLsf)
      end

      obj.iBlts       = iBlts;
      obj.ssid        = ssid;
      obj.isAchg      = isAchg;
      obj.iCalibTimeL = iCalibTimeL;
      obj.iCalibTimeH = iCalibTimeH;
      obj.iLsf        = iLsf;
    end



  end



end
