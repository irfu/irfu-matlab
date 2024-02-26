%
% Default algorithm for selecting output dataset filename for BPCIs. Meant to be
% used as an argument to bicas.tools.batch.autocreate_one_SWM_BPCI(),
% unless the caller wants to use their own customized function.
%
% NOTE: Not obvious whether to use CDAG, which version.
% NOTE: A wrapper function can set settings by adding them to the end of the
% argument list. Migt not be able to do this using an anonymous function though.
%
%
% ARGUMENTS
% =========
% BpciInputDsmdArray
%       Column array of solo.adm.DSMD. DSMDs for BPCI input datasets. These are
%       used for setting the time the output dataset filename.
% varargin
%       Arguments as interpreted by
%       irf.utils.interpret_settings_args().
% --
% Determined by argument "createOutputFilenameFunc" to function
% bicas.tools.batch.autocreate_one_SWM_BPCI().
%
%
% RETURN VALUES
% =============
% filename
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-15.
%
function filename = default_get_BPCI_output_filename(...
  outputDsi, BpciInputDsmdArray, ...
  versionStr, varargin)

% PROPOSAL: Automatic test code.
%
% PROPOSAL: INPUT_DSI_FOR_OUTPUT_TIME --> argument
%   PROPOSAL: Only use first matching DSI in INPUT_DSI_FOR_OUTPUT_TIME.
%       PRO: More general. Less constraint on SWMs.
%
% PROPOSAL: Add argument for parent directory.
%   CON: Setting parent directory is a different task. Different from
%        defining filenaming conventions.
%   PRO: Common use case.

% List of DSIs. Only matching input datasets are used for deriving
% time range in output filename. May contain unused DSIs.
INPUT_DSI_FOR_OUTPUT_TIME = { ...
  'SOLO_L1_RPW-LFR-SBM1-CWF', ...
  'SOLO_L1_RPW-LFR-SBM2-CWF', ...
  'SOLO_L1_RPW-LFR-SURV-CWF', ...
  'SOLO_L1_RPW-LFR-SURV-SWF', ...
  'SOLO_L1_RPW-TDS-LFM-CWF', ...
  'SOLO_L1_RPW-TDS-LFM-RSWF', ...
  ...
  'SOLO_L1R_RPW-LFR-SBM1-CWF-E', ...
  'SOLO_L1R_RPW-LFR-SBM2-CWF-E', ...
  'SOLO_L1R_RPW-LFR-SURV-CWF-E', ...
  'SOLO_L1R_RPW-LFR-SURV-SWF-E', ...
  'SOLO_L1R_RPW-TDS-LFM-CWF-E', ...
  'SOLO_L1R_RPW-TDS-LFM-RSWF-E', ...
  ...
  'SOLO_L2_RPW-LFR-SURV-CWF-E'};

DEFAULT_SETTINGS.unoffBasenameExtension = [];       % No extension.
DEFAULT_SETTINGS.isCdagPolicy           = false;    % true/false/<other>

Settings = irf.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
irf.assert.struct(Settings, fieldnames(Settings), {})



% ASSERTIONS
assert(ischar(outputDsi))
assert(isa(BpciInputDsmdArray, 'solo.adm.DSMD'))
%assert(iscell(cohbCa))
assert(ischar(versionStr))
%assert(numel(BpciInputDsmdArray) == numel(cohbCa))
assert(islogical(Settings.isCdagPolicy) || isnumeric(Settings.isCdagPolicy), ...
  'Illegal Settings.isCdagPolicy.')



% Identify exactly one BPCI INPUT DSMD which shall be used for determining
% time interval for OUTPUT dataset.
[~, iDsmd] = intersect({BpciInputDsmdArray.datasetId}, INPUT_DSI_FOR_OUTPUT_TIME);
assert(isscalar(iDsmd), ...
  'Can not determine exactly one input DSI to use for determining output filename time interval.')
InputDsmd = BpciInputDsmdArray(iDsmd);



%==========================================
% Set variables in output dataset filename
%==========================================
R = struct();
R.isCdag         = logical(Settings.isCdagPolicy);
R.datasetId      = outputDsi;
R.versionStr     = versionStr;
R.dsicdagCase    = 'lower';
R.unoffExtension = Settings.unoffBasenameExtension;

% Set date vector(s), depending on time range, effectively selecting filename
% format for the dataset.
dt1 = InputDsmd.dt1;
dt2 = InputDsmd.dt2;
if dt2 <= (dt1 + caldays(1))
  % CASE: (dt1,dt2) covers less than or equal to a calendar day.
  % ==> Use filenaming format yymmdd (no begin-end; just the calendar day).

  % IMPLEMENTATION NOTE: Using center of day so that can handle DSMD time
  % boundaries that are slightly on the wrong side of midnight, e.g. if
  % deriving DSMD from file content.
  %   Ex: LFR-SURV-CWF test file
  %       solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200213_V07.cdf begins
  %       at Epoch = 2020-02-12T23:59:53.305345024
  dtMiddle = dt1 + (dt2-dt1)/2;
  R.dateVec  = [dtMiddle.Year, dtMiddle.Month, dtMiddle.Day];
else
  % CASE: (dt1,dt2) covers more than one day.
  % ==> use filenaming format yymmddThhmmss-yymmddThhmmss.
  R.dateVec1 = datevec(dt1);
  R.dateVec2 = datevec(dt2);
end

%================================
% Create output dataset filename
%================================
filename = solo.adm.create_dataset_filename(R);
end
