%
% YK 2020-10-15: Do not use any unofficial basename extension for the IRF
% pipeline. This is necessary for irfu-matlab's automatic zVar & dataset
% finding.
%
%
% ARGUMENTS
% =========
% BpciInputDsmdArray
%       DSMD column array for INPUT datasets for the (one) relevant BPCI.
% PreexistingOutputLvDsmdArray
%       DSMD column array containing only the latest versions (LV)
%       of preexisting datasets which shall be used for determining output
%       dataset version number.
% fnVerAlgorithm
%       String constant. Represents selected filename (FN) version algorithm.
%
%
% RETURN VALUES
% =============
% filePath
%       Path to output file, including filename.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function filePath = get_BPCI_output_path2(...
  BpciInputDsmdArray, PreexistingOutputLvDsmdArray, ...
  outputDsi, fnVerAlgorithm, outputDir, outputIsCdag)

% PROPOSAL: Test for day with leap second.
%
% PROPOSAL: Rename PreexistingOutputLvDsmdArray
%   dataset version
%   consideration
%   preexisting
%   basis (for)
%   competitors
%
% PROPOSAL: Not set output directory. The caller should do that.
%
% PROBLEM: SBM1/2 has L1R datasets on a different filenaming (time) format
%   Ex:
%     solo_L1R_rpw-lfr-sbm1-cwf-e-cdag_20240201T025448-20240201T030848_V02.cdf
%     ==>
%     solo_L2_rpw-lfr-sbm1-cwf-e_20240201_V01.cdf                         (incorrect)
%     solo_L2_rpw-lfr-sbm1-cwf-e_20240201T025448-20240201T030848_V01.cdf  (correct)
%   PROBLEM: DSMDs do not store the filenaming format.
%   PROPOSAL: Separate list of input DSIs which should yield this format.
%   PROPOSAL: Use length of time interval. -- IMPLEMENTED
%
% PROPOSAL: Separate function for converting a selected reference input
%           dataset filename into output dataset filename. -- IMPLEMENTED

assert(isa(BpciInputDsmdArray,           'solo.adm.DSMD') && iscolumn(BpciInputDsmdArray))
assert(isa(PreexistingOutputLvDsmdArray, 'solo.adm.DSMD') && iscolumn(PreexistingOutputLvDsmdArray))
assert(ischar(outputDsi))
assert(ischar(fnVerAlgorithm))
assert(islogical(outputIsCdag))

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



% Identify exactly one BPCI INPUT DSMD which shall be used for determining
% time interval for OUTPUT dataset.
[~, iRefDsmd] = intersect({BpciInputDsmdArray.datasetId}, INPUT_DSI_FOR_OUTPUT_TIME);
assert(isscalar(iRefDsmd), ...
  'Can not determine exactly one input DSI to use for determining output filename time interval.')
InputRefDsmd = BpciInputDsmdArray(iRefDsmd);

dt1 = InputRefDsmd.dt1;
dt2 = InputRefDsmd.dt2;

versionNbr = get_output_version(...
  dt1, dt2, ...
  outputDsi, fnVerAlgorithm, PreexistingOutputLvDsmdArray);

fileName = get_BPCI_output_filename2(...
  dt1, dt2, ...
  outputDsi, outputIsCdag, versionNbr);

filePath = fullfile(outputDir, fileName);
end



function versionNbr = get_output_version(...
  dt1, dt2, outputDsi, fnVerAlgorithm, PreexistingOutputLvDsmdArray)

%=================================================================
% Identify highest version number of pre-existing output datasets
%=================================================================
PreexistingOutputLvDsmdArray = solo.adm.filter_DSMD_DATASET_ID(...
  PreexistingOutputLvDsmdArray, {outputDsi});

if ~isempty(PreexistingOutputLvDsmdArray)
  % NOTE: Command only works for non-empty array.
  bKeep1 = ([PreexistingOutputLvDsmdArray.dt1] == dt1);
  bKeep2 = ([PreexistingOutputLvDsmdArray.dt2] == dt2);
  bKeep = bKeep1 & bKeep2;
else
  bKeep = zeros(0, 1);
end
PreexistingOutputLvDsmd = PreexistingOutputLvDsmdArray(bKeep);   % Empty or scalar.

if isscalar(PreexistingOutputLvDsmd)
  versionNbr = PreexistingOutputLvDsmd.versionNbr;
elseif isempty(PreexistingOutputLvDsmd)
  versionNbr = NaN;
else
  error('Could not find zero or one matching DSMD in PreexistingOutputLvDsmdArray.')
end

%============================================
% Determine version number of output dataset
%============================================
switch(fnVerAlgorithm)

  case 'ABOVE_HIGHEST_USED'
    if isnan(versionNbr)
      versionNbr = 1;
    else
      versionNbr = versionNbr + 1;
    end

  case 'HIGHEST_USED'
    if isnan(versionNbr)
      versionNbr = 1;
    else
      versionNbr = versionNbr;   % Nonsense command for clarity.
    end

  otherwise
    error('Illegal fnVerAlgorithm="%s".', fnVerAlgorithm)
end
end



function outputFileName = get_BPCI_output_filename2(...
  dt1, dt2, outputDsi, outputIsCdag, versionNbr)

UNOFF_BASENAME_EXTENSION = [];

%==========================================
% Set variables in output dataset filename
%==========================================
R = [];
R.datasetId          = outputDsi;
R.versionStr         = sprintf('%02i', versionNbr);
R.isCdag             = outputIsCdag;
R.dsicdagCase        = 'lower';
R.unoffExtension     = UNOFF_BASENAME_EXTENSION;
R.dateVec1           = datevec(dt1);
R.dateVec2           = datevec(dt2);

% Set date vector, depending on time range, effectively selecting filename
% time interval format for the dataset.

% if dt2 <= (dt1 + caldays(1))
if is_midnight(dt1) && is_midnight(dt2) && (dt2 == dt1 + caldays(1))
  % CASE: (dt1,dt2) covers exactly one calender day.
  % ==> Use filenaming format yymmdd (no begin-end; just the calendar day).
  %
  % NOTE/BUG: This is technically somewhat imperfect since input dataset
  % would
  % solo_L1R_rpw-lfr-sbm1-cwf-e-cdag_20240201T000000-20240202T000000_V02.cdf
  % would yield an output dataset filename on the YYYYMMDD format. This
  % should be unlikely.

  R.timeIntervalFormat = 'DAY';
else
  % CASE: (dt1,dt2) does not exactly cover one day one day.
  % ==> use filenaming format yymmddThhmmss-yymmddThhmmss.
  R.timeIntervalFormat = 'SECOND_TO_SECOND';
end



%================================
% Create output dataset filename
%================================
outputFileName = solo.adm.dsfn.create_dataset_filename(R);
end



function isMidnight = is_midnight(dt)
isMidnight = dateshift(dt, 'start', 'day') == dt;
end
