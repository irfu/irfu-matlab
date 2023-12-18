%
% Class that represents how to convert one particular QRCID into modifications
% of quality ZVs.
%
% NOTE: Class does not include the QRCID itself.
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
    %   QRC/QRCID
    %       CON: Class refers to translation of QRCID to quality variables modifications
    %   quality ZVs
    %   quality ZVs modification
    %   QrcSetting
    %
    properties(SetAccess=immutable)
        % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
        %
        % NOTE: This is interpretation is in compliance with how the ZV
        % QUALITY_FLAG is supposed to be set/updated.
        QUALITY_FLAG

        % Bits that should be set in either ZV "L2_QUALITY_BITMASK" or
        % "L3_QUALITY_BITMASK". The context in which this class is used
        % determines which.
        Lx_QUALITY_BITMASK
    end



    methods(Access=public)



        function obj = QrcSetting(QUALITY_FLAG, Lx_QUALITY_BITMASK)

            assert(bicas.utils.validate_ZV_QUALITY_FLAG(QUALITY_FLAG))
            assert(isa(Lx_QUALITY_BITMASK, 'uint16'))

            obj.QUALITY_FLAG       = QUALITY_FLAG;
            obj.Lx_QUALITY_BITMASK = Lx_QUALITY_BITMASK;
        end



    end



end
