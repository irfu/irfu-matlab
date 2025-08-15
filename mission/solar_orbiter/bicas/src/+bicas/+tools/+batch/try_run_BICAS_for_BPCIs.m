%
% Try run BICAS for array of specified BPCIs. If the input datasets can not be
% found immediately before a call to BICAS, then that BPCI is skipped. This is
% to prevent the code from confusing a "genuine" BICAS error (which the code can
% not interpret) with BICAS merely not being able to read missing input
% datasets.
%
%
% RETURN VALUE
% ============
% BpcsArray
%       Column array of BPCSs for the BICAS calls made.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function BpcsArray = try_run_BICAS_for_BPCIs(...
  Bpa, BpciArray, configFile, automountTriggerPathsCa, bicasSettingsArgsCa)

% PROPOSAL: Better name. Omit "try".
%   CON: "try" implies handling failure.
%       CON: Does imply catching exception, only that exception may occur.

assert(isa(Bpa,       'bicas.tools.batch.BicasProcessingAccessAbstract'))
assert(isa(BpciArray, 'bicas.tools.batch.BicasProcessingCallInfo'))
assert(iscolumn(BpciArray))
assert(iscell(bicasSettingsArgsCa))

% BpcsArray = bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1);
nBpci  = numel(BpciArray);
BpcsCa = cell(nBpci, 1);

%################################
% Loop over BPCIs/calls to BICAS
%################################
% NOTE: Loop can be run both using (1) "for", and (2) "parfor" (i.e.
% parallelized).

% NOTE: Using "parfor" seems to work (and output CDFs are identical), except
% that:
% (1) The stdout from bicas.tools.batch contains a mix of stdout from different
%     BICAS instances being executed in parallel.
% (2) The stdout from bicas.tools.batch contains some seemingly irrelevant
%     warnings (1) from MATLAB.
% --
% NOTE: A quick test (irony, network-mounted input
% /data/solo/remote/data/L2/lfr_wf_e/2025/06/*cwf*_2025060*) reduced the
% execution time used 197 s --> 121 s. Shorter test showed no improvement,
% probably due to the processing pool startup time.

for iBpci = 1:nBpci
% parfor iBpci = 1:nBpci
  Bpci = BpciArray(iBpci);

  Bpcs = execute_BICAS_using_BPCI(...
    Bpa, Bpci, configFile, bicasSettingsArgsCa, automountTriggerPathsCa, iBpci);

  if isempty(Bpcs)
    % CASE: BICAS has not been run (called)
    % -------------------------------------
    % NOTE: There shall be no BPCS when BICAS has not run. BPCS stores any
    %       BICAS error, but not that it has not run at all.
    BpcsCa{iBpci, 1} = bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1);
  else
    BpcsCa{iBpci, 1} = Bpcs;
  end
end

%------------------------------------------------------------------------------
% Convert BpcsCa (some elements contain []) --> BpcsArray (no elements are [])
%------------------------------------------------------------------------------
% IMPLEMENTATION NOTE: Doing this in a separate for loop to potentially make it
% possible to use parfor for the previous for loop.
% IMPLEMENTATION NOTE: Can not use cell2mat() since it does not support cell
% arrays of objects.
BpcsArray = bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1);
for iBpcs = 1:nBpci
  Bpcs = BpcsCa{iBpcs};

  if ~isempty(Bpcs)
    BpcsArray(end+1, 1) = Bpcs;
  end
end
end    % try_run_BICAS_for_BPCIs()



function Bpcs = execute_BICAS_using_BPCI(...
  Bpa, Bpci, configFile, bicasSettingsArgsCa, automountTriggerPathsCa, iBpci)

assert(numel(Bpci.outputsArray) >= 1)

argsCa = bicas.tools.batch.BPCI_to_BICAS_call_args(Bpci);

% NOTE: Uses first output file (i.e. no other output files) to name the
%       log file.
outputPath   = Bpci.outputsArray(1).path;
timestampStr = datestr(now, 'YYYY-mm-ddTHH.MM.SS');
logFile      = sprintf('%s.%s.log', outputPath, timestampStr);

% IMPLEMENTATION NOTE: Overriding BICAS settings wrt. SWMs used for
% identifying can NOT be done here only since the exact set of SWMs is needed
% earlier for grouping input datasets.
argsCa(end+1:end+2) = {'--config',     configFile};
argsCa(end+1:end+2) = {'--log-matlab', logFile};
argsCa              = [argsCa(:); bicasSettingsArgsCa(:)];

% Trigger automounts
% ==================
% IMPLEMENTATION NOTE: The BPCIs contain canonical paths due to the use of
% dir() when generating paths to all files under a specified path using
% bicas.tools.batch.get_file_paths(). Can therefore not trigger automounts
% using the paths used in the call to BICAS directly (i.e. in the BPCI).
irf.fs.trigger_automounts(automountTriggerPathsCa)

%========================================
% Check if input datasets are available,
% otherwise skip and try the next BPCI
%========================================
[inputPathsValid, firstInvalidPath] = input_dataset_paths_valid(Bpci);
if ~inputPathsValid
  warning([...
    'Can not find previously identified input file needed for a call', ...
    ' to BICAS. This may e.g. be due to that', ...
    ' (1) it has been moved during the execution of this script, or', ...
    ' (2) the disk was temporarily inaccessible (NAS automount problem).', ...
    ' It will be searched for again in the next pass.\n', ...
    '    firstInvalidPath="%s"'], ...
    firstInvalidPath)
  Bpcs = [];
  return
end

% IMPLEMENTATION NOTE: It is useful to have log messages specifying when the
% call to BICAS begins and ends to identify overlapping BICAS instances in
% stdout when using parallelization (parfor).
fprintf('######################\n')
fprintf('CALL BICAS (iBpci=%i)\n', iBpci)
fprintf('######################\n')
fprintf(...
  'CALLING BICAS WITH THE FOLLOWING ARGUMENTS: %s\n', ...
  ['"', strjoin(argsCa, sprintf('"\n    "')), '"'])

%############
% CALL BICAS
%############
errorCode = Bpa.bicas_main(argsCa{:});

fprintf('################################\n')
fprintf('RETURNING FROM BICAS (iBpci=%i)\n', iBpci)
fprintf('################################\n')

% Create BPCS for the completed call.
Bpcs = bicas.tools.batch.BicasProcessingCallSummary(Bpci, errorCode);
end    % execute_BICAS_using_BPCI()



% Check whether the input datasets for one particular BPCI actually exist. Meant
% to be used just before launching BICAS with that BPCI to prevent launching
% BICAS for input datasets that have just been moved (in particular moved to
% subdirectory former_versions/ due to syncing ROC data while executing this
% code).
%
% NOTE: Does NOT work on an BPCI array since the check has to be done JUST
% BEFORE launching BICAS.
%
%
% RETURN VALUES
% =============
% valid
%       Logical. Whether Bpci input datasets can be found or not.
% firstInvalidPath
%       String. Path. Only set when valid==true.
%
function [valid, firstInvalidPath] = input_dataset_paths_valid(Bpci)
assert(isscalar(Bpci))

valid            = true;    % Value until set to the opposite.
firstInvalidPath = [];
for i = 1:numel(Bpci.inputsArray)
  path = Bpci.inputsArray(i).path;
  if ~exist(path, 'file')
    firstInvalidPath = path;
    valid            = false;
  end
end
end    % input_dataset_paths_valid()
