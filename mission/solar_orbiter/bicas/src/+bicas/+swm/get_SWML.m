%
% Get official SWML used by BICAS (as opposed to any future test SWML).
%
%
% NOTE
% ====
% NOTE: There is no SWM for generating VHT datasets(!). VHT processing is
%       therefore not represented here.
% NOTE: The SWMs implemented in the body of BICAS must always be compatible
%       with the SWMs specified here.
%
%
% RETURN VALUES
% =============
% Official bicas.swm.SoftwareModeList object (~singleton) that is actually used
% by BICAS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function Swml = get_SWML(Bso)
% PROPOSAL: Re-implement (top-level) hard-coded constants by setting
%   multiple redundant 1D(?) vectors that covers every case. Then
%   set various cases by assigning constants to many elements using
%   MATLAB syntax. One index representing: Combination of
%   DSI+Skeleton_Version (both pipelines, LFR+TDS, HK+SCI),
%   every element contains data for that dataset. Must use
%   combination DSI+Skeleton_Version to potentially cover old
%   versions. Manipulate and set multiple elements smoothly by using
%   vectors for indices.
%   Ex: Vectors to set: skeletonVersionVector, SBMx_SURV_vector,
%       CWF_SWF_vector, output dataset level, vectors for
%       human-readable description string(s) (e.g. modeStr)
%   Ex: Vectors with indices for dataset in/for: either pipeline,
%       science or HK, LFR or TDS, latest versions or
%       backward-compatibility versions.
%   --
%   TODO-DEC: Above describes input & output (?) data sets.
%       How relates to s/w modes?
%
% PROPOSAL: Merge LFR and TDS loops.
% PROPOSAL: Split constructor into separate functions, one per s/w
%           mode (or at least the big ones)
%   L1/L1R --> L2 mode
%   L2     --> L3
%
% PROPOSAL: Replace BSO with arguments for each used setting.
%   CON: BSO (and L) is used by production functions. Can not
%        at all eliminate easily.

% Arrays with constants.
% {1} = S/w modes (science) L1R-->L2
% {2} = S/w modes (science) L1 -->L2
inputDatasetLevelList     = {'L1R'};
inputDashEList            = {'-E'};
swmSuffixList             = {''};
swmPurposeAmendmList      = {''};
if Bso.get_fv('SWM.L1-L2_ENABLED')
  inputDatasetLevelList{end+1} = 'L1';
  inputDashEList{end+1}        = '';
  swmSuffixList{end+1}         = '_L1';
  swmPurposeAmendmList{end+1}  = ' UNOFFICIAL wrt. ROC.';
end



% Define function which interprets (replaces) specific substrings.
% "strmod" = string modify, "g"=global
strmodg = @(s, iInputLevel) bicas.utils.strrep_many(s, ...
  '<InLvl>',              inputDatasetLevelList{iInputLevel}, ...
  '<I-E>',                inputDashEList{iInputLevel}, ...
  '<SWM suffix>',         swmSuffixList{iInputLevel}, ...
  '<SWM purpose amendm>', swmPurposeAmendmList{iInputLevel});

% Input definitions that are reused multiple times.
HK_INPUT_DEF = bicas.swm.InputDataset(...
  'in_hk', 'SOLO_HK_RPW-BIA', 'HK_cdf');
CUR_INPUT_DEF = bicas.swm.InputDataset(...
  'in_cur', 'SOLO_L1_RPW-BIA-CURRENT', 'CUR_cdf');



% Define arrays of data used for generating s/w modes definitions.
% One component per pair of s/w modes L1/L1R-->L2
LFR_SWM_DATA = struct(...
  'SBMx_SURV',       {'SBM1', 'SBM2', 'SURV', 'SURV'}, ...
  'CWF_SWF',         {'CWF',  'CWF',  'CWF',  'SWF' }, ...
  'modeStr',         {...
  'selective burst mode 1', ...
  'selective burst mode 2', ...
  'survey mode', ...
  'survey mode' ...
  }, ...
  'outputSkeletonVersion', {'16', '16', '16', '16'});
TDS_SWM_DATA = struct(...
  'CWF_RSWF',              {'CWF', 'RSWF'}, ...
  'outputSkeletonVersion', {'16',  '16'});



SwmList = bicas.swm.SoftwareMode.empty(0, 1);
% Iterate over [L1R] (one component), or [L1, L1R]...
for iInputLevel = 1:numel(inputDatasetLevelList)

  %==============================================
  % Iterate over the "fundamental" LFR S/W modes
  %==============================================
  for iSwm = 1:length(LFR_SWM_DATA)

    % Define local string modification function.
    strmod = @(s) bicas.utils.strrep_many(strmodg(s, iInputLevel), ...
      '<SBMx/SURV>',  LFR_SWM_DATA(iSwm).SBMx_SURV, ...
      '<C/SWF>',      LFR_SWM_DATA(iSwm).CWF_SWF, ...
      '<mode str>',   LFR_SWM_DATA(iSwm).modeStr);

    SciInputDataset = bicas.swm.InputDataset(...
      'in_sci', ...
      strmod('SOLO_<InLvl>_RPW-LFR-<SBMx/SURV>-<C/SWF><I-E>'), ...
      'SCI_cdf');

    SciOutputDataset = bicas.swm.OutputDataset(...
      'out_sci', ...
      strmod('SOLO_L2_RPW-LFR-<SBMx/SURV>-<C/SWF>-E'), ...
      'SCI_cdf', ...
      strmod('LFR L2 <C/SWF> science electric <mode str> data'), ...
      strmod(...
      ['RPW LFR L2 <C/SWF> science electric', ...
      ' (potential difference) data in <mode str>,', ...
      ' time-tagged']), ...
      LFR_SWM_DATA(iSwm).outputSkeletonVersion);

    SwmList(end+1) = bicas.swm.SoftwareMode(...
      bicas.proc.L1L2.LfrSwmProcessing(...
      SciInputDataset.dsi, SciOutputDataset.dsi), ...
      strmod('LFR-<SBMx/SURV>-<C/SWF>-E<SWM suffix>'), ...
      strmod(...
      ['Generate <SBMx/SURV> <C/SWF> electric field', ...
      ' L2 data (potential difference)', ...
      ' from LFR <InLvl> data.<SWM purpose amendm>']), ...
      [SciInputDataset, CUR_INPUT_DEF, HK_INPUT_DEF], ...
      [SciOutputDataset]);
  end

  %==============================================
  % Iterate over the "fundamental" TDS S/W modes
  %==============================================
  for iSwm = 1:numel(TDS_SWM_DATA)

    if strcmp(TDS_SWM_DATA(iSwm).CWF_RSWF, 'RSWF') ...
        && strcmp(inputDatasetLevelList{iInputLevel}, 'L1')
      % CASE: TDS RSWF
      % Exclude SWM since BICAS can not currently (2023-10-09) read
      % TDS RSWF L1 datasets.!
      continue
    end

    % Define local string modification function.
    strmod = @(s) strrep(strmodg(s, iInputLevel), ...
      '<C/RSWF>', TDS_SWM_DATA(iSwm).CWF_RSWF);

    SciInputDataset = bicas.swm.InputDataset(...
      'in_sci', ...
      strmod('SOLO_<InLvl>_RPW-TDS-LFM-<C/RSWF><I-E>'), ...
      'SCI_cdf');

    SciOutputDataset = bicas.swm.OutputDataset(...
      'out_sci', ...
      strmod('SOLO_L2_RPW-TDS-LFM-<C/RSWF>-E'), ...
      'SCI_cdf', ...
      strmod('LFR L2 <C/RSWF> science electric LF mode data'), ...
      strmod(...
      ['RPW TDS L2 <C/RSWF> science electric (potential', ...
      ' difference) data in LF mode, time-tagged']), ...
      TDS_SWM_DATA(iSwm).outputSkeletonVersion);

    SwmList(end+1) = bicas.swm.SoftwareMode(...
      bicas.proc.L1L2.TdsSwmProcessing(...
      SciInputDataset.dsi, SciOutputDataset.dsi), ...
      strmod('TDS-LFM-<C/RSWF>-E<SWM suffix>'), ...
      strmod(...
      ['Generate <C/RSWF> electric field L2 data', ...
      ' (potential difference)', ...
      ' from TDS LF mode <InLvl> data.<SWM purpose amendm>']), ...
      [SciInputDataset, CUR_INPUT_DEF, HK_INPUT_DEF], ...
      [SciOutputDataset]);
  end
end    % for iInputLevel = 1:numel(inputDatasetLevelList)



if Bso.get_fv('SWM.L2-L2_CWF-DSR_ENABLED')
  SciInputDataset = bicas.swm.InputDataset(...
    'in_sci', ...
    'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
    'OSR_cdf');

  SciOutputDataset = bicas.swm.OutputDataset(...
    'out_dsr', ...
    'SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND', ...
    'DSR_cdf', ...
    'LFR L2 CWF science electric survey mode data, downsampled', ...
    ['RPW LFR L2 CWF science electric (potential difference)', ...
    ' data in survey mode, time-tagged, downsampled'], ...
    '02');

  % NOTE: Function handle: Argument rctDir is not used, but is
  % needed for the interface.
  SwmList(end+1) = bicas.swm.SoftwareMode(...
    bicas.proc.L2L2.LfrDsrSwmProcessing(), ...
    'LFR-SURV-CWF-E-DSR', ...
    ['Generate downsampled version of LFR L2 SURV CWF', ...
    ' science electric (potential difference) data.', ...
    ' NOTE: This is an unofficial s/w mode.'], ...
    [SciInputDataset], ...
    [SciOutputDataset]);
end



if Bso.get_fv('SWM.L2-L3_ENABLED')
  SciInputDataset = bicas.swm.InputDataset(...
    'in_sci', ...
    'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
    'LFR-SURV-CWF-E_cdf');

  EfieldOutputDataset = bicas.swm.OutputDataset(...
    'out_efield', ...
    'SOLO_L3_RPW-BIA-EFIELD', ...
    'EFIELD_OSR_cdf', ...
    'BIAS L3 science electric field vector data', ...
    'RPW BIAS L3 science electric field vector data, time-tagged', ...
    '05');

  ScpotOutputDataset = bicas.swm.OutputDataset(...
    'out_scpot', ...
    'SOLO_L3_RPW-BIA-SCPOT', ...
    'SCPOT_OSR_cdf', ...
    'BIAS L3 science spacecraft potential data', ...
    'RPW BIAS L3 science spacecraft potential data, time-tagged', ...
    '05');

  DensityOutputDataset = bicas.swm.OutputDataset(...
    'out_density', ...
    'SOLO_L3_RPW-BIA-DENSITY', ...
    'DENSITY_OSR_cdf', ...
    'BIAS L3 science plasma density data', ...
    'RPW BIAS L3 science plasma density data, time-tagged', ...
    '06');

  EfieldDsrOutputDataset = bicas.swm.OutputDataset(...
    'out_efield_dsr', ...
    'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
    'EFIELD_DSR_cdf', ...
    'BIAS L3 downsampled science electric field vector data', ...
    ['RPW BIAS L3 downsampled science electric', ...
    ' field vector data, time-tagged'], ...
    '05');

  ScpotDsrOutputDataset = bicas.swm.OutputDataset(...
    'out_scpot_dsr', ...
    'SOLO_L3_RPW-BIA-SCPOT-10-SECONDS', ...
    'SCPOT_DSR_cdf', ...
    'BIAS L3 downsampled science spacecraft potential data', ...
    ['RPW BIAS L3 downsampled science spacecraft', ...
    ' potential data, time-tagged'], ...
    '05');

  DensityDsrOutputDataset = bicas.swm.OutputDataset(...
    'out_density_dsr', ...
    'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS', ...
    'DENSITY_DSR_cdf', ...
    'BIAS L3 downsampled science plasma density data', ...
    ['RPW BIAS L3 downsampled science plasma', ...
    ' density data, time-tagged'], ...
    '06');

  % NOTE: Function handle: Arguments rctDir, NsoTable are not
  % used, but are needed for the interface.
  SwmList(end+1) = bicas.swm.SoftwareMode(...
    bicas.proc.L2L3.L3OsrDsrSwmProcessing, ...
    'BIA-EFIELD-SCPOT-DENSITY', ...
    ['Generate L3 electric field vector, spacecraft', ...
    ' potential, and density data', ...
    ' incl. additional downsampled versions.', ...
    ' NOTE: This is an unofficial s/w mode.'], ...
    [SciInputDataset], ...
    [EfieldOutputDataset, ...
    EfieldDsrOutputDataset, ...
    ScpotOutputDataset, ...
    ScpotDsrOutputDataset, ...
    DensityOutputDataset, ...
    DensityDsrOutputDataset]);
end

Swml = bicas.swm.SoftwareModeList(SwmList);

end
