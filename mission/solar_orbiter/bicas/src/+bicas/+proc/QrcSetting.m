%
% Class that represents how to convert one particular QRCID into modifications
% of quality ZVs.
%
% NOTE: The class does not include the QRCID itself.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSetting
  % PROPOSAL: Better name
  %   setting
  %       CON: Could be confused with more proper settings such as the NSO
  %            table itself. The information stored in the object is more like
  %            a constant.
  %   interpretation
  %   configuration
  %   behaviour
  %   entry
  %   policy
  %   entry
  %   action
  %   QRC, QRCID
  %   quality ZVs
  %   quality ZVs modification
  %   --
  %   QRCC = QRC Configuration
  %   QRCS = QRC Setting -- IMPLEMENTED
  %   --
  %   PROPOSAL: Abbreviation
  %
  % PROPOSAL: Use only keyword arguments with default values.
  %   PRO: Common to not set all values.
  %     Ex: Hardcoding tests.
  %   CON: Longer calls.



  properties(SetAccess=immutable)
    % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
    %
    % NOTE: This is interpretation is in compliance with how the ZV
    % QUALITY_FLAG is supposed to be set/updated.
    L2_QUALITY_FLAG
    L3_QUALITY_FLAG_density

    % Bits (bitmask) that should be set in either ZV "L2_QUALITY_BITMASK" or
    % "L3_QUALITY_BITMASK".
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    L2_QUALITY_BITMASK
    L3_QUALITY_BITMASK_density

    % Column array (set) of unique SSIDs for which L2 voltage samples should be
    % blanked.
    l2VoltageFvSsidAr

    % Column array (set) of unique antenna numbers (1-3) for which L2 current
    % samples should be blanked.
    l2CurrentFvIantAr
  end



  methods(Access=public)



    function obj = QrcSetting(A)
      arguments
        A.L2_QUALITY_FLAG            uint8  = bicas.const.QUALITY_FLAG_MAX
        A.L3_QUALITY_FLAG_density    uint8  = bicas.const.QUALITY_FLAG_MAX
        A.L2_QUALITY_BITMASK         uint16 = uint16(0)
        A.L3_QUALITY_BITMASK_density uint16 = uint16(0)
        A.l2VoltageFvSsidAr          uint8  = uint16.empty(0, 1)
        A.l2CurrentFvIantAr                 = zeros(0, 1)
      end

      assert(bicas.utils.validate_ZV_QUALITY_FLAG(A.L2_QUALITY_FLAG))
      assert(bicas.utils.validate_ZV_QUALITY_FLAG(A.L3_QUALITY_FLAG_density))
      obj.L2_QUALITY_FLAG         =               A.L2_QUALITY_FLAG;
      obj.L3_QUALITY_FLAG_density =               A.L3_QUALITY_FLAG_density;

      assert(isa(                       A.L2_QUALITY_BITMASK,         'uint16'))
      assert(isa(                       A.L3_QUALITY_BITMASK_density, 'uint16'))
      obj.L2_QUALITY_BITMASK          = A.L2_QUALITY_BITMASK;
      obj.L3_QUALITY_BITMASK_density  = A.L3_QUALITY_BITMASK_density;

      assert(iscolumn(                     A.l2VoltageFvSsidAr))
      assert(bicas.proc.L1L2.const.is_SSID(A.l2VoltageFvSsidAr))
      irf.assert.number_set(               A.l2VoltageFvSsidAr)
      obj.l2VoltageFvSsidAr  =             A.l2VoltageFvSsidAr;

      assert(iscolumn(         A.l2CurrentFvIantAr))
      assert(all(ismember(     A.l2CurrentFvIantAr, [1, 2, 3])))
      irf.assert.number_set(   A.l2CurrentFvIantAr)
      obj.l2CurrentFvIantAr  = A.l2CurrentFvIantAr;
    end



  end



end
