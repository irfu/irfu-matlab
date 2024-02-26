%
% Convert a BPCI into a sequence of arguments that can be used for calling
% BICAS.
%
%
% Initially created 2020-03-02 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function argsCa = BPCI_to_BICAS_call_args(Bpci)
% PROPOSAL: Refactor into BPCI method.
% PROPOSAL: No hardcoded CLI_OPTION_PREFIX
%   PROPOSAL: CLI_OPTION_PREFIX as argument.
%   PROPOSAL: CLI_OPTION_PREFIX using BICAS constant.

    CLI_OPTION_PREFIX = '--';

    assert(isa(Bpci, 'bicas.tools.batch.BicasProcessingCallInfo'))

    argsCa = {Bpci.swmCliOption};

    for i = 1:numel(Bpci.inputsArray)
        In            = Bpci.inputsArray(i);
        argsCa{end+1} = [CLI_OPTION_PREFIX, In.cohb];
        argsCa{end+1} =                     In.path;
    end

    for i = 1:numel(Bpci.outputsArray)
        Out           = Bpci.outputsArray(i);
        argsCa{end+1} = [CLI_OPTION_PREFIX, Out.cohb];
        argsCa{end+1} = Out.path;
    end
end
