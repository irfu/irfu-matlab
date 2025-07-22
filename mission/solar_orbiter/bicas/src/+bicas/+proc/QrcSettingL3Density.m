%
% QRCS for producing L3 density datasets.
%
% NOTE: The class does not include the QRCID itself.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL3Density < bicas.proc.QrcSetting



  properties(SetAccess=immutable)
    % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
    % NOTE: This is interpretation is in compliance with how the ZV
    % QUALITY_FLAG is supposed to be set/updated.
    QUALITY_FLAG

    % Bits (bitmask) that should be set in either ZV "L2_QUALITY_BITMASK" or
    % "L3_QUALITY_BITMASK" (density).
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    L3_QUALITY_BITMASK
  end



  methods(Access=public)



    function obj = QrcSettingL3Density(A)
      arguments
        A.QUALITY_FLAG       uint8  = bicas.const.QUALITY_FLAG_MAX
        A.L3_QUALITY_BITMASK uint16 = uint16(0)
      end

      assert(bicas.utils.validate_ZV_QUALITY_FLAG(A.QUALITY_FLAG))
      obj.QUALITY_FLAG =                          A.QUALITY_FLAG;

      assert(isa(               A.L3_QUALITY_BITMASK, 'uint16'))
      obj.L3_QUALITY_BITMASK  = A.L3_QUALITY_BITMASK;
    end



  end



end
