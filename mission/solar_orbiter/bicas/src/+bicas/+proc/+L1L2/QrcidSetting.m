%
% Class that represents how to interpret a particular QRCID.
%
% NOTE: Does not include the QRCID itself.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcidSetting
    % PROPOSAL: Better name
    %   setting
    %       CON: Could be confused with more proper settings such as the NSO
    %            table itself. The information stored in the object is more like
    %            a constant.
    %   interpretation
    %   behaviour
    %   entry
    %   policy
    %   entry
    %   action



    properties(SetAccess=immutable)
        QUALITY_FLAG
        % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
        %
        % NOTE: This is interpretation is in compliance with how the ZV
        % QUALITY_FLAG is supposed to be set.

        L2_QUALITY_BITMASK
        % Bits that should be set in ZV "L2_QUALITY_BITMASK".
    end



    methods(Access=public)

        function obj = QrcidSetting(QUALITY_FLAG, L2_QUALITY_BITMASK)
            assert(bicas.utils.validate_ZV_QUALITY_FLAG(QUALITY_FLAG))
            assert(isa(L2_QUALITY_BITMASK, 'uint16'))

            obj.QUALITY_FLAG       = QUALITY_FLAG;
            obj.L2_QUALITY_BITMASK = L2_QUALITY_BITMASK;
        end

    end



end
