%
% Class that represents how to convert one particular QRC (QRCID into
% modifications of quality ZVs, blanking of data etc. for a particular DSI.
% Implemented by subclasses.
%
% NOTE: The user of the subclasses should know their interface, and therefore
% this abstract superclass does not *NOT* define any interface (sic!). This is
% since different information is needed for different forms of processing.
% Ex: ZVs L2_QUALITY_BITMASK vs L3_QUALITY_BITMASK.
% Ex: ZV QUALITY_FLAG is technically shared between all output datasets, but
%     could be regarded as qualititatively different for different forms of
%     processing(i.e. warrants separate variables)
%
% NOTE: The class does not itself the QRCID or DSI themselves.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) QrcSetting



  methods(Access=public)



    function obj = QrcSetting()
    end



  end



end
