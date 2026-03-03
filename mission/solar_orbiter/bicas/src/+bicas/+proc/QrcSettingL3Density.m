%
% QRCS for producing L3 *DENSITY* datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingL3Density < bicas.proc.QrcSetting



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Bits (bitmask) that should be set in ZV "L3_QUALITY_BITMASK" (density).
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    l3qbm
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingL3Density(A)
      arguments
        A.qfl       = bicas.const.qrc.QFL_MAX
        A.l3qbm     = bicas.const.qrc.LxQBM_NONE
        A.gaCaveats = string.empty(0, 1);
      end

      obj@bicas.proc.QrcSetting(qfl=A.qfl, gaCaveats=A.gaCaveats);

      assert(isa( A.l3qbm, 'uint16'))
      obj.l3qbm = A.l3qbm;
    end



  end



end
