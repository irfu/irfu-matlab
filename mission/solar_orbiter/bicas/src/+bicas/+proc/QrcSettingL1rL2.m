%
% Extension of superclass, specific for L1/L1R-->L2 processing.
%
% NOTE: Does not contain L2_QUALITY_BITMASK since that is still contained in
% the superclass (sic!) since it is referenced in generic code for both
% L1/L1R-->L2 and L2-->L3 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL1rL2 < bicas.proc.QrcSetting



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Array (set) of SSIDs for which voltage samples should be blanked.
    voltageSamplesFvSsidAr
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingL1rL2(...
        QUALITY_FLAG, L2_QUALITY_BITMASK, voltageSamplesFvSsidAr)
      obj = obj@bicas.proc.QrcSetting(QUALITY_FLAG, L2_QUALITY_BITMASK);

      assert(iscolumn(voltageSamplesFvSsidAr))
      assert(bicas.proc.L1L2.const.is_SSID(voltageSamplesFvSsidAr))
      assert(numel(unique(voltageSamplesFvSsidAr)) == numel(voltageSamplesFvSsidAr))

      obj.voltageSamplesFvSsidAr = voltageSamplesFvSsidAr;
    end



  end    % methods(Access=public)



end
