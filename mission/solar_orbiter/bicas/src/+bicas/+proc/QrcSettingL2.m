%
% QRCS for producing official (OSR) L2 datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL2 < bicas.proc.QrcSetting



  properties(SetAccess=immutable)
    % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
    % NOTE: This is interpretation is in compliance with how the ZV
    % QUALITY_FLAG is supposed to be set/updated.
    QUALITY_FLAG

    % Bits (bitmask) that should be set in ZV "L2_QUALITY_BITMASK" or
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    L2_QUALITY_BITMASK

    % Column array (set) of unique SSIDs for which L2 voltage samples should be
    % blanked.
    voltageFvSsidAr

    % Column array (set) of unique antenna numbers (1-3) for which L2 current
    % samples should be blanked.
    currentFvIantAr
  end



  methods(Access=public)



    function obj = QrcSettingL2(A)
      arguments
        A.QUALITY_FLAG       uint8  = bicas.const.QUALITY_FLAG_MAX
        A.L2_QUALITY_BITMASK uint16 = uint16(0)
        A.voltageFvSsidAr    uint8  = uint16.empty(0, 1)
        A.currentFvIantAr           = zeros(0, 1)
      end

      assert(bicas.utils.validate_ZV_QUALITY_FLAG(A.QUALITY_FLAG))
      obj.QUALITY_FLAG            =               A.QUALITY_FLAG;

      assert(isa(                   A.L2_QUALITY_BITMASK, 'uint16'))
      obj.L2_QUALITY_BITMASK      = A.L2_QUALITY_BITMASK;

      assert(iscolumn(                     A.voltageFvSsidAr))
      assert(bicas.proc.L1L2.const.is_SSID(A.voltageFvSsidAr))
      irf.assert.number_set(               A.voltageFvSsidAr)
      obj.voltageFvSsidAr  =             A.voltageFvSsidAr;

      assert(iscolumn(         A.currentFvIantAr))
      assert(all(ismember(     A.currentFvIantAr, [1, 2, 3])))
      irf.assert.number_set(   A.currentFvIantAr)
      obj.currentFvIantAr  = A.currentFvIantAr;
    end



  end



end
