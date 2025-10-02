%
% QRCS for producing official (OSR) L2 datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL2 < bicas.proc.QrcSetting



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Bits (bitmask) that should be set in ZV "L2_QUALITY_BITMASK" or
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    l2qbm

    % Column array (set) of unique SSIDs for which L2 voltage samples should be
    % blanked, *BEFORE* reconstruction.
    voltageFvSsidAr

    % Column array (set) of unique antenna numbers (1-3) for which L2 current
    % samples should be blanked.
    currentFvIantAr
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingL2(A)
      arguments
        A.qfl             = bicas.const.qrc.QFL_MAX
        A.l2qbm           = bicas.const.qrc.LxQBM_NONE
        A.voltageFvSsidAr = uint8.empty(0, 1)
        A.currentFvIantAr = zeros(0, 1)
        A.gaCaveats       = string.empty(0, 1);
      end

      obj@bicas.proc.QrcSetting(qfl=A.qfl, gaCaveats=A.gaCaveats);

      assert(isa( A.l2qbm, 'uint16'))
      obj.l2qbm = A.l2qbm;

      assert(iscolumn(                     A.voltageFvSsidAr))
      assert(bicas.proc.L1L2.const.is_SSID(A.voltageFvSsidAr))
      irf.assert.number_set(               A.voltageFvSsidAr)
      obj.voltageFvSsidAr  =               A.voltageFvSsidAr;

      assert(iscolumn(      A.currentFvIantAr))
      assert(all(ismember(  A.currentFvIantAr, [1, 2, 3])))
      irf.assert.number_set(A.currentFvIantAr)
      obj.currentFvIantAr = A.currentFvIantAr;
    end



  end



end
