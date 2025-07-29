%
% Class that represents how to convert one particular QRC (QRCID into
% modifications of quality ZVs, blanking of data etc. for an implicit group of
% DSIs. Implemented by subclasses which represent different groups of DSIs.
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
  % PROPOSAL: Include QUALITY_FLAG.


  methods(Access=public)



    function obj = QrcSetting()
    end



  end



end
