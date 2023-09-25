%
% Immutable class which instances represent the destination of a signal, i.e.
% either:
% (1) an ASR (to determine how to store the signal in the dataset), or
% (2) "nowhere" (when mux mode is unknown).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalDestinationId



    properties(GetAccess=public, Constant)
        C = bicas.proc.L1L2.SignalDestinationId.init_const();
    end



    properties(SetAccess=immutable, GetAccess=public)
        value
    end



    methods(Access=public)



        % Constructor
        function obj = SignalDestinationId(value)
            assert(isequal(value, 'Nowhere') || isa(value, 'bicas.proc.L1L2.AntennaSignalId'))

            obj.value = value;
        end



    end    % methods(Access=public)



    methods(Access=private, Static)



        function C = init_const()
            C = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants( ...
                @(asid) (bicas.proc.L1L2.SignalDestinationId(asid)));

            C.NOWHERE = bicas.proc.L1L2.SignalDestinationId('Nowhere');
        end



    end



end
