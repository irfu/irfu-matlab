%
% Autocreate all BPCIs which can be derived from specified paths.
%
% NOTE: If the reference directory is simultaneously specified as an input path
% (optional), then that is included here too.
%
%
% ARGUMENTS
% =========
% get_BPCI_output_path_fh
%       Function handle.
%       path = @(outputDsi, BpciInputDsmdArray, cohbCa).
%       Determines file paths for BPCI output datasets.
%
%
% RETURN VALUES
% =============
% BpciArray
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function BpciArray = autocreate_input_BPCIs(...
  InputDsmdArray, get_BPCI_output_path_fh, SwmArray, currentDatasetExtensionDays)

% IMPLEMENTATION NOTE: Setting will probably eventually be abolished by the
% function solo.adm.group_sort_DSMD_versions() that uses it, but it is kept
% here as a global constant (instead of an argument) for clarity.
SETTINGS_sortWrtFormerVersionsDir = false;



%=====================================
% Find all DSMDs in INPUT directories
%=====================================
t = tic();
fprintf('INPUT paths: #Datasets, all versions (DSMDs): %i\n', numel(InputDsmdArray))
fprintf('SPEED: bicas.tools.batch.autocreate_input_BPCIs(): Time to obtain input DSMDs: %.1f [s] wall time\n', toc(t));



%============================
% Filter out latest versions
%============================

% Must do BEFORE filtering latest version.
InputDsmdArray = solo.adm.convert_DSMD_DATASET_ID_to_SOLO(InputDsmdArray);

% NOTE: Might be slow.
t = tic();
InputDsmdArray = solo.adm.group_sort_DSMD_versions(...
  InputDsmdArray, 'latest', ...
  'sortWrtFormerVersionsDir', SETTINGS_sortWrtFormerVersionsDir);
fprintf('SPEED: bicas.tools.batch.autocreate_input_BPCIs(): solo.adm.group_sort_DSMD_versions(): %.1f [s] wall time\n', toc(t));

% NOTE: Must do AFTER filtering latest version. The function extends the
% DSMD time intervals which would make solo.adm.group_sort_DSMD_versions()
% fail to group the datasets correctly.
InputDsmdArray = solo.adm.extend_last_CURRENT_DSMD(...
  InputDsmdArray, currentDatasetExtensionDays);
fprintf('INPUT paths: #Datasets (DSMDs), latest versions: %i\n', numel(InputDsmdArray))



%=============================================
% Find all BPCIs, based on available datasets
%=============================================
t = tic();
BpciArray = bicas.tools.batch.autocreate_many_BPCIs(...
  InputDsmdArray, SwmArray, ...
  get_BPCI_output_path_fh);
fprintf('SPEED: bicas.tools.batch.autocreate_input_BPCIs(): autocreate_many_BPCIs(): %.1f [s] wall time\n', toc(t));
fprintf('#BPCIs, found total: %3i\n', numel(BpciArray))
end
