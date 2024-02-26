%
% Each instance summarizes a completed BICAS call, including error code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BicasProcessingCallSummary



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        Bpci

        % BICAS error code
        errorCode
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = BicasProcessingCallSummary(Bpci, errorCode)

            assert(isa(Bpci, 'bicas.tools.batch.BicasProcessingCallInfo'))
            assert(isnumeric(errorCode))

            obj.Bpci      = Bpci;
            obj.errorCode = errorCode;
        end



    end    % methods(Access=public)



end
