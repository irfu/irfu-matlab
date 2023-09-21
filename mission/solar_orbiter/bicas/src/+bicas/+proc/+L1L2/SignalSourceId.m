%
% Immutable class which instances represent the source of a signal, i.e.
% either:
% (1) an ASR, or
% (2) various special cases.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalSourceId
% PROPOSAL: Properties for "category" and "antennas" that mirror the
%           bicas.proc.L1L2.AntennaSignalId stored inside the object in most
%           cases. Error if other.



    properties(GetAccess=public, Constant)
        C = bicas.proc.L1L2.SignalSourceId.init_const()
    end



    properties(SetAccess=immutable, GetAccess=public)
        value
    end



    methods(Access=public)

        % Constructor
        function obj = SignalSourceId(value)
            assert(...
                isequal(value, '2.5V Ref') || ...
                isequal(value, 'GND'     ) || ...
                isequal(value, 'Unknown' ) || ...
                isa(value, 'bicas.proc.L1L2.AntennaSignalId') ...
            )

            obj.value = value;
        end

        function isAsr = is_ASR(obj)
            isAsr = isa(obj.value, 'bicas.proc.L1L2.AntennaSignalId');
        end

    end    % methods(Access=public)



    methods(Access=private, Static)

        function SSID = init_const()
            SSID = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants( ...
                @(asid) (bicas.proc.L1L2.SignalSourceId(asid)));

            SSID.REF25V   = bicas.proc.L1L2.SignalSourceId('2.5V Ref');
            SSID.GND      = bicas.proc.L1L2.SignalSourceId('GND');
            SSID.UNKNOWN  = bicas.proc.L1L2.SignalSourceId('Unknown');
        end

    end



end
