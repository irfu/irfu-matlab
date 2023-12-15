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
    %
    % PROPOSAL: Replace QUALITY_FLAG and *_QUALITY_BITMASK with containers.Map:
    %           key--> QUALITY_FLAG and *_QUALITY_BITMASK value.
    %   PRO: Can simultaneously specify different behaviour different types of output datasets.
    %   TODO-DEC: What should map key be?
    %       PROPOSAL: DSI
    %           CON: Too many identical entries.
    %       PROPOSAL: Informal group of output datasets.
    %           Ex: L2, L3_density
    %   PROPOSAL: Class method with DSI as argument returns only relevant
    %             values.
    %       CON: Can not handle (hypothetically) qualitatively different types
    %            of quality variables for different groups of output datasets,
    %            if method output should always be on the same format.
    %           Ex: Density has L3_QUALITY_BITMASK but EFIELD and SCPOT do not
    %               (I think).
    %           Ex: Hypothetical: Varying number of bits in *_QUALITY_BITMASK in
    %               different output datasets.
    %           CON: There is no (real) such case.



    properties(SetAccess=immutable)
        % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
        %
        % NOTE: This is interpretation is in compliance with how the ZV
        % QUALITY_FLAG is supposed to be set.
        QUALITY_FLAG

        % Bits that should be set in ZV "L2_QUALITY_BITMASK".
        L2_QUALITY_BITMASK
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
