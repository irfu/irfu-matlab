%
% Run multiple passes over input datasets, making it potentially possible to
% process the input of previous passes.
%
% NOTE: See bicas.tools.batch.main() for detailed documentation on
% the exact "algorithm" etc.
%
%
% ARGUMENTS
% =========
% SwmArray
%       SWM array.
%       NOTE: Use for identifying input datasets and setting output datasets
%       etc.
%
%
% RETURN VALUES
% =============
% BpcsAllArray
%       Array of BPCSs. Metadata on all calls made to BICAS, in cronological
%       order.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [BpcsAllArray] = run_BICAS_all_passes(...
  Bpa, bicasSettingsArgsCa, configFile, ...
  outputDir, referenceDir, ...
  inputPathsCa, fnVerAlgorithm, outputIsCdag, ...
  SwmArray, Settings)

% NOTE: Algorithm documentation in bicas.tools.batch.main().
%
% PROPOSAL: Better name.
%
% PROPOSAL: Convert iteration to function. ~run_BICAS_one_pass().
%   PRO: Can convert get_BPCI_output_path_fh to regular function.
%   PRO: Size of code can be allowed to be larger.
%
% PROPOSAL: Never convert ref.dir. --> DSDMs.
%   PRO: Only filenames are used.
%   CON: Must filter out the datasets from other files though.
%       CON: Not really since they are only used for matching against output
%            dataset filenames.
%           CON: Inefficient.
%
% PROPOSAL: tpdFilenamesCa should contain paths, not filenames.
%   PRO: DSMDs store paths.
%   CON: TPD datasets always refer to output dir. paths which need to
%        compared with (for set ops. diff+union) ref. dir. paths. Those
%        operations must be done on filenames anyway.
%   PROPOSAL: DSMDs for TPD using only filenames.
%
% PROPOSAL: Log (print) input files, output files.
%   CON: Could be many.
%       CON-PROPOSAL: Optional for debugging.
%       CON: Useful for logging in general.
%           CON: BICAS own log messages include files in & out.

DEBUG_ENABLED = false;
% DEBUG_ENABLED = true;

assert(isa(Bpa, 'bicas.tools.batch.BicasProcessingAccessAbstract'))
assert(iscell(inputPathsCa) && iscolumn(inputPathsCa), 'inputPathsCa is not a column cell array.')

% List of output datasets which BICAS has TRIED to create, successfully or
% not. If BICAS failed to produce them, then the corresponding files do not
% exist, but the code should still not try to generate them again since
% BICAS will likely just fail again and for every future attempt forever.
tpdFilenamesCa = cell(0,1);
% BPCSs for all passes so far.
BpcsAllArray   = bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1);



iPass = 1;
while true

  fprintf('########################\n')
  fprintf('BEGIN BICAS BATCH PASS %i\n', iPass)
  fprintf('########################\n')

  %=====================
  % Get reference DSMDs
  %=====================
  if ~isempty(referenceDir)
    [filePathsCa,   ~] = bicas.tools.batch.get_file_paths({referenceDir});
    [RefDsmdArray,  ~] = solo.adm.paths_to_DSMD_array(filePathsCa);
  else
    RefDsmdArray = solo.adm.DSMD.empty(0, 1);
  end
  if DEBUG_ENABLED
    log_DSMD_array('Reference directory (only recognized datasets)', RefDsmdArray)
  end

  %=====================================================================
  % Function handle to function which creates OUTPUT DATASET PATHS,
  % including filenames, given relevant BPCI, datasets to consider when
  % determining output data versions etc.
  %=====================================================================
  % IMPLEMENTATION NOTE: Removing TPD file names to make sure that any
  % previously generated file would be re-generated by re-running (what is
  % effectively) the same BPCI.
  preexistingOutputFilenamesCa = setdiff(...
    bicas.tools.batch.DSMDs_to_filenames(RefDsmdArray), ...
    tpdFilenamesCa);

  % NEW ALGORITHM: FASTER
  % ---------------------
  % IMPLEMENTATION NOTE: Converts filenames to DSMDs to be able to
  % efficiently extract only the latest versions (LV). Should possibly
  % convert preexistingOutputFilenamesCa to DSMDs or paths to make
  % this natural.
  PreexistingOutputDsmdArray = solo.adm.paths_to_DSMD_array(preexistingOutputFilenamesCa(:));
  PreexistingOutputLvDsmdArray = solo.adm.group_sort_DSMD_versions(...
    PreexistingOutputDsmdArray, ...
    'latest', 'sortWrtFormerVersionsDir', false);

  get_BPCI_output_path_fh = ...
    @(outputDsi, BpciInputDsmdArray) ( ...
    bicas.tools.batch.get_BPCI_output_path2(...
    BpciInputDsmdArray, PreexistingOutputLvDsmdArray, ...
    outputDsi, fnVerAlgorithm, ...
    outputDir, outputIsCdag));

  %====================================================
  % Autocreate all possible BPCIs based on input paths
  %====================================================
  if isempty(referenceDir)
    referenceDirCa = cell(0, 1);
  else
    referenceDirCa = {referenceDir};
  end
  [filePathsCa,    ~] = bicas.tools.batch.get_file_paths(inputPathsCa);
  [InputDsmdArray, ~] = solo.adm.paths_to_DSMD_array(filePathsCa);
  if DEBUG_ENABLED
    log_DSMD_array('Input directories (only recognized datasets)', InputDsmdArray)
  end

  BpciInputArray = bicas.tools.batch.autocreate_input_BPCIs(...
    InputDsmdArray, get_BPCI_output_path_fh, SwmArray, ...
    Settings.currentDatasetExtensionDays);
  fprintf('Number of possible BPCIs:         %3i\n', numel(BpciInputArray))

  %=======================================================================
  % Find out which subset of BPCIs that should actually be run
  % ----------------------------------------------------------
  % NOTE: Takes reference directory into account but behaviour depends on
  %       fnVerAlgorithm in get_BPCI_output_path_fh.
  %   fnVerAlgorithm = 'HIGHEST_USED':
  %       Output dataset filenames have same version as highest
  %       counterpart in ref. dir., if there is one. ==> Filename
  %       collision. ==> Excluded
  %   fnVerAlgorithm = 'ABOVE_HIGHEST_USED':
  %       Output dataset filenames have a higher version than highest
  %       version counterpart in ref. dir., if there is one.
  %       ==> Never filename collision. ==> Included/kept.
  %=======================================================================
  doNotNeedToGenerateFilenamesCa = union(...
    bicas.tools.batch.DSMDs_to_filenames(RefDsmdArray), ...
    tpdFilenamesCa);
  BpciRunArray = bicas.tools.batch.filter_BPCIs_to_run(...
    BpciInputArray, doNotNeedToGenerateFilenamesCa);
  fprintf('Number of BPCIs that will be run: %3i\n', numel(BpciRunArray))

  if isempty(BpciRunArray)
    % CASE: NO MORE BPCIs TO RUN
    break
  end

  %=============================================================================
  % Try run BICAS for selected BPCIs
  % --------------------------------
  % Skip BPCIs for which not all input datasets exist just before
  % execution.
  %
  % NOTE: automountTriggerPathsCa includes inputPathsCa which can be very long
  % if it refers to explicit datasets. Could possibly slow down execution.
  % PROPOSAL: Only include the first N paths of inputPathsCa in automountTriggerPathsCa.
  %   CON: It is not yet known to be a problem.
  %=============================================================================
  automountTriggerPathsCa = [{configFile; outputDir; referenceDir}; inputPathsCa];
  BpcsPassArray = bicas.tools.batch.try_run_BICAS_for_BPCIs(...
    Bpa, BpciRunArray, configFile, automountTriggerPathsCa, bicasSettingsArgsCa);

  BpcsAllArray = [BpcsAllArray; BpcsPassArray];



  % BpcsAllArray --> Output dataset filenames = tpdFilenamesCa
  % (= TPD filenames for all passes so far.)
  tpdFilenamesCa = get_BPCSs_output_filenames(BpcsAllArray);

  iPass = iPass + 1;
end



% Assign return values.
becArray = [BpcsAllArray.errorCode];
nTpd     = numel(tpdFilenamesCa);
end



% Function for clarifying the code.
function filenamesCa = get_BPCSs_output_filenames(BpcsArray)
filenamesCaCa = arrayfun(...
  @(Bpci) (Bpci.get_output_filenames()), ...
  [BpcsArray(:).Bpci]', ...
  'UniformOutput', false);

filenamesCa = cat(1, filenamesCaCa{:});
end



% Log DSMD array by listing paths, adding header/title and empty rows before and
% after.
%
% Can be used for debugging, if not for logging permanently.
function log_DSMD_array(titleStr, DsmdArray)
fprintf('\n')
fprintf('\n%s\n', titleStr)
fprintf('%s\n', repmat('=', 1, numel(titleStr)))

nDsmd = numel(DsmdArray);
for i = 1:nDsmd
  Dsmd = DsmdArray(i);
  fprintf('    %s\n', Dsmd.path)
end
% IMPLEMENTATION NOTE: If there are zero DSMDs, then it is really important
% to print something instead of zero rows of paths since that is confusing
% (looks like a potential bug to the user).
fprintf('Total: %i\n', nDsmd)

fprintf('\n')
end
