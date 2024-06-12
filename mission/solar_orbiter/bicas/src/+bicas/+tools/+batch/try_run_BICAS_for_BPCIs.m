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

BpcsArray = bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1);

for iBpci = 1:numel(BpciArray)
  Bpci = BpciArray(iBpci);

  assert(numel(Bpci.outputsArray) >= 1)

  argsCa = bicas.tools.batch.BPCI_to_BICAS_call_args(Bpci);

  % NOTE: Uses first output file (i.e. no other output files) to name the
  %       log file.
  outputPath   = Bpci.outputsArray(1).path;
  timestampStr = datestr(now, 'YYYY-mm-ddTHH.MM.SS');
  logFile      = sprintf('%s.%s.log', outputPath, timestampStr);

  % IMPLEMENTATION NOTE: Overriding BICAS settings wrt. SWMs used for
  % identifying can NOT be done here only since the exact set of SWMs
  % is needed earlier for grouping input datasets.
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
    continue
  end

  %====================================
  % Run BICAS
  % Create BPCS for the completed call
  %====================================
  fprintf(...
    'CALLING BICAS WITH THE FOLLOWING ARGUMENTS: %s\n', ...
    ['"', strjoin(argsCa, sprintf('"\n    "')), '"'])

  %############
  % CALL BICAS
  %############
  errorCode = Bpa.bicas_main(argsCa{:});

  Bpcs = bicas.tools.batch.BicasProcessingCallSummary(...
    Bpci, errorCode);

  BpcsArray(end+1, 1) = Bpcs;
end
end



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
end
