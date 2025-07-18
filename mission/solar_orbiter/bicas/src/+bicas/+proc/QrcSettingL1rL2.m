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
  % PROPOSAL: Use only keyword arguments with default values.
  %   PRO: Common to not set all values.
  %     Ex: Hardcoding tests.
  %   CON: Longer calls.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Column array (set) of unique SSIDs for which voltage samples should be blanked.
    voltageSamplesFvSsidAr

    % Column array (set) of unique antenna numbers (1-3) for which current samples should be
    % blanked.
    currentSamplesFvIAntAr
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingL1rL2(...
        QUALITY_FLAG, L2_QUALITY_BITMASK, ...
        voltageSamplesFvSsidAr, ...
        currentSamplesFvIantAr)

      obj = obj@bicas.proc.QrcSetting(QUALITY_FLAG, L2_QUALITY_BITMASK);

      assert(iscolumn(voltageSamplesFvSsidAr))
      assert(bicas.proc.L1L2.const.is_SSID(voltageSamplesFvSsidAr))
      assert(numel(unique(voltageSamplesFvSsidAr)) == numel(voltageSamplesFvSsidAr))
      obj.voltageSamplesFvSsidAr = voltageSamplesFvSsidAr;

      assert(iscolumn(currentSamplesFvIantAr))
      assert(all(ismember(currentSamplesFvIantAr, [1, 2, 3])))
      assert(numel(unique(currentSamplesFvIantAr)) == numel(currentSamplesFvIantAr))
      obj.currentSamplesFvIAntAr = currentSamplesFvIantAr;
    end



  end    % methods(Access=public)



end
