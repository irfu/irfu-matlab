%
% Keep only those BPCIs which should actually be run. Always keep
% BPCI if any output file is missing in both reference directory (unless empty)
% and tpdFilenameCa.
%
% NOTE: Uses exact output filenames in argument to determine which BPCIs should
% not be run. Therefore output filenames have to match for
% (1) CDAG/non-CDAG,
% (2) version number
% (3) unofficial basename suffix, if used.
%
%
% ARGUMENTS
% =========
% doNotNeedToGenerateFilenamesCa
%       Cell array of dataset filenames for datasets.
%
%
% RETURN VALUE
% ============
% BpciArray
%       Input argument BpciArray but only that subset of elements where at least
%       one of the output datasets has the same filename as on in
%       generateIfMissingFilenamesCa.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function BpciArray = filter_BPCIs_to_run(BpciArray, doNotNeedToGenerateFilenamesCa)
    % PROPOSAL: Better name.
    %   Something more generic with filtering.
    %   filter, keep, remove
    %   output datasets
    %   BPCI

    % ASSERTIONS
    assert(isa(BpciArray, 'bicas.tools.batch.BicasProcessingCallInfo'))
    assert(iscolumn(BpciArray))
    assert(iscell(doNotNeedToGenerateFilenamesCa))

    %============================================================
    % Filter BpciArray:
    % Only keep BPCIs for which at least one datasets is missing
    %============================================================
    bKeep = false(size(BpciArray));
    for iBpci = 1:numel(BpciArray)
        outputFilenameCa = BpciArray(iBpci).get_output_filenames();

        % Keep BPCI if at least one of its output datasets is missing.
        bKeep(iBpci) = ~all(ismember(outputFilenameCa, doNotNeedToGenerateFilenamesCa));
    end

    BpciArray = BpciArray(bKeep, 1);
end
