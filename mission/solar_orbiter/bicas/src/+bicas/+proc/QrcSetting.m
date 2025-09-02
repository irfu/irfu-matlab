%
% Class that represents how to convert one particular QRC (QRCID) into
% modifications of quality ZVs, blanking of data etc. for an implicit group of
% DSIs. Implemented by subclasses which represent different groups of output
% DSIs. This implicitly means that the same QRC might apply to multiple
% instances (of different subclasses), e.g. for both L2 and L3 output.
%
% NOTE: The user of the subclasses should know their interface, and therefore
% this abstract superclass does not *NOT* define any interface (sic!). The
% reason is that different information is needed for different forms of
% processing.
% Ex: ZVs L2_QUALITY_BITMASK vs L3_QUALITY_BITMASK.
% Ex: ZV QUALITY_FLAG is technically shared between all output datasets, but
%     could be regarded as qualititatively different for different forms of
%     processing (i.e. warrants separate variables)
%
% NOTE: This superclass is still referenced by QRCSM by for assertions.
% NOTE: The class does not itself the QRCID or DSI themselves.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) QrcSetting
  % PROPOSAL: Include QUALITY_FLAG/qfl.
  %   NOTE: Currently not a property in bicas.proc.QrcSettingL3.   /2025-09-02



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % GA "CAVEATS". Column array of strings.
    gaCaveats
  end



  methods(Access=public)



    function obj = QrcSetting(A)
      arguments
        A.gaCaveats = string.empty(0, 1);
      end

      assert(iscolumn(A.gaCaveats) & isstring(A.gaCaveats))
      obj.gaCaveats = A.gaCaveats;
    end



  end



end
