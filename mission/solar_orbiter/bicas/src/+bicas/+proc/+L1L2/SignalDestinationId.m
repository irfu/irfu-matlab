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
        % ASID object or empty.
        Asid
        
        % Whether destination is "nowhere", i.e. the signal does not have a
        % destination and should be ignored.
        isNowhere
    end



    methods(Access=public)



        % Constructor
        function obj = SignalDestinationId(value)
            if isa(value, 'bicas.proc.L1L2.AntennaSignalId')
                obj.Asid      = value;
                obj.isNowhere = false;
            elseif isequal(value, 'Nowhere')
                obj.Asid      = [];
                obj.isNowhere = true;
            else
                error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
            end
        end



    end    % methods(Access=public)



    methods(Access=private, Static)



        function C = init_const()
            C = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants( ...
                @(Asid) (bicas.proc.L1L2.SignalDestinationId(Asid)));

            C.NOWHERE = bicas.proc.L1L2.SignalDestinationId('Nowhere');
        end



    end



end
