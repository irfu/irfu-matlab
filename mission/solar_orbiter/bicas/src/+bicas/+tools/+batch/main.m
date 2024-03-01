%
% Code for batch processing of datasets with BICAS.
%
% NOTE: There is a separate wrapper script for calling this function from a bash
% script.
%
%
% NOTES
% =====
% IMPORTANT NOTE: Function intentionally continues on BICAS error (on-zero error
% code). Will then not try the same BICAS call again. Will only call BICAS if
% all input datasets exist immediately prior to the call.
% --
% NOTE: Will try to only use latest versions of input datasets.
% NOTE: Output dataset filenames are determined by
%       bicas.tools.batch.default_get_BPCI_output_filename().
% NOTE: Will produce one log file per BICAS run.
% NOTE: Will itself not configure settings to enable extra SWMs. Use the BICAS
%       config file to enable/disable SWMs as far as BICAS permits it.
%
% NOTE: Relies on BICAS code for
%   (1) processing
%   (2) obtaining SWM definitions (and potentially uses BICAS SETTINGS to obtain
%       this not yet implemented).
%       NOTE: Therefore effectively relies on an unofficial BICAS interface!
% --
% IMPLEMENTATION NOTE: If one wants to specify
%   ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH,
%   ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH,
% then one needs a system-customized BICAS config file. Therefore, the
% implementation REQUIRES the caller to specify a config file.
% --
% NOTE: The time needed to identify DSMDs and BPCIs can be slow if applied to
% many datasets. This problem will grow as the mission collects more and more
% data. Therefore logs time consumption for these parts in anticipation of havng
% to speed this up.
% --
% NOTE: If one uses the output directory as a source directory, then the
% datasets in the output directory will also be used as input for BICAS in
% subsequent rounds of dataset searching+processing. This makes it possible to
% process L1R-->L2-->L3 with a single call.
% --
% YK 2020-10-15: Do not use any unofficial basename extension for the IRF
% pipeline. This is necessary for irfu-matlab's automatic zVar & dataset
% finding.
%
%
% ALGORITHM FOR WHICH DATASETS THAT WILL BE PRODUCED
% ==================================================
% DEFINITION: ONE PASS
%   Let BICAS produce all datasets which can be produced from the datasets in
%   the given input paths (given all existing BICAS SWMs), as long as there is
%   no corresponding output dataset (any version) in the reference directory.
%   Different combinations of referenceDir and modeStr therefore determine which
%   datasets are produced and which versions the output datasets will have.
%   NOTE: The implementation uses filename collisions between output datasets
%   and reference directory datasets in a way which is not what one would expect
%   from the description above, but the resulting behaviour is the same.
% --
% Run pass after pass in a loop, until no new datasets can be produced.
% NOTE: Incrementing dataset versions are always done relative to reference
% directory.
%
%
% ALGORITHM FOR HOW TO IDENTIFY BPCIs
% ===================================
% See bicas.tools.batch.autocreate_many_BPCIs().
% NOTE: Algorithm requires input datasets to not overlap in time for each
% DSI separately. Can therefore not handle SOLO_L2_RPW-LFR-SBM1-CWF-E
%
%
% HOW TO PROCESS: EXAMPLES
% ========================
% PROCESS EVERYTHING (potentially overwrite datasets)
%   referenceDir = []
%   modeStr      = V01 (or NEW_VERSIONS)
% PROCESS EVERYTHING, INCREMENT VERSIONS
%   referenceDir = directory with older datasets (can be output directory)
%   modeStr      = NEW_VERSION
% ONLY PROCESS NON-PREEXISTING DATASETS
%   referenceDir = directory with older datasets (can be output directory)
%   modeStr      = V01
% PROCESS L1/L1R+HK --> L2 --> L3 IN ONE CALL
%   referenceDir = directory with older datasets (can be output directory)
%   modeStr      = V01
%   inputPathsCa : Includes referenceDir    ## IMPORTANT
%
%
% ARGUMENTS
% =========
% bicasConfigFile
%       Path to BICAS' config file.
% outputIsCdag
%       Logical. Whether output datasets should be CDAG.
% modeStr : String constant.
%       'V01'
%           Only generate datasets not already represented in
%           referenceDir (any version). Output datasets always have version V01.
%       'NEW_VERSION'
%           Output datasets are set to having the first unused version number
%           over the highest pre-existing ones among the corresponding output
%           datasets in the reference directory.
%           NOTE: Primarily intended to be used for L3 deliveries to ROC.
% outputDir
%       Directory to create new datasets in.
% referenceDir
%       Directory path. Directory used for determining which datasets have
%       already been produced. Used by argument "modeStr". Will be searched
%       recursively.
%       --
%       SPECIAL CASE: Empty <=> Empty reference directory.
% inputPathsCa : Cell array
%       Arbitrary number of paths (arguments) to individual datasets and/or
%       directories. Datasets are searched for (recursively) under directories.
%       Other files are ignored. The found datasets are used for finding
%       combinations of datasets which can be used to generate datasets, based
%       on filenaming conventions.
% varargin
%       Settings as interpreted by
%       irf.utils.interpret_settings_args().
%
%
% INTERNAL NAMING CONVENTIONS
% ===========================
% BEC  : BICAS Error Code
% HTPD : Have or Tried Produce Dataset(s).
% TPD  : Tried Produce Dataset(s). Output dataset which the code actually tried
%        to produce by calling BICAS (successfully or not), as opposed to output
%        datasets that could not be produced due to removed input datasets
%        detected before calling BICAS.
% --
% NOTE: Code also uses a small number of abbreviations defined in BICAS.
%
%
% Initially created 2020-03-02 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function main(...
  bicasConfigFile, outputIsCdag, modeStr, ...
  outputDir, referenceDir, inputPathsCa, varargin)
% PROPOSAL: Better name
%   ~bicas
%       CON: Clashes with package name
%   ~batch
%       CON: Clashes with package name
%           CON: Is sort of in line with Python's convention of main module
%                in a package sharing name with the (leaf-most) package name.
%
% PROPOSAL: More natural to use bicas.Logger, bicas.utils.Timekeeper.
%
% PROPOSAL: Be able to optionally limit to specified SWM(s).
%   PROPOSAL: Use settings.
%       CON: How handle from bash?
%   TODO-DEC: How specify set of SWMs?
%       PROPOSAL: Specify set of SWMs to permit (include)
%       PROPOSAL: Specify set of SWMs ignore (exclude)
%       TODO-DEC: Permit non-existent SWMs?
%           NOTE: The SWMs which BICAS sees depends on BICAS settings.
%           PROPOSAL: Not permit
%               PRO: More rigorous. Less change of misspelling.
%               CON: Can not write wrappers, aliases with permanent lists of
%                    SWMs.
%
% PROPOSAL: Filter DSIs (only keep those needed for SWMs) before sorting
%           by version and before searching for BPCIs.
%   PRO: Likely faster.
%
% PROPOSAL: Rename modeStr --> outputVersionModeStr
%   CON: modeStr does not only determine the output version algorithm. V01
%        only generates datasets not already in reference directory
%        whereas NEW_VERSION always generates new datasets.
%
% PROPOSAL: Redefine algorithm.
%   NOTE: Technically, the reference directory is used for two different
%       things:
%       (1) Determine output dataset versions.
%       (2) Determine whether to output datasets.
%           NOTE: modeStr=NEW_VERSION will always always produce output
%           dataset candidates which never collide with reference datasets.
%   --
%   PROPOSAL: Argument/flag for whether to permit overwriting datasets.
%       PROPOSAL: Not argument visible in outside interface. Could be a
%                 function of "mode" for at least the short term.
%       TODO-DEC: Would need policy on how to handle SWMs/BPCIs which
%           simultaneously output multiple files, but for which there is a
%           collision only for some.
%           NOTE: Always write all files or no file. + Overwrite disallowed.
%                 ==> Any collision implies writing no file.
%       PRO: Clarity.
%       PRO: Flexibility.
%           CON: More to test.
%   PROPOSAL:
%       Generate BPCIs with output dataset versions using reference
%       directory + versioning algorithm (as now).
%       Run every BPCI which does not only produce datasets pre-existing in
%       reference directory.
%
% NOTE/~BUG/~PROBLEM: Algorithm requires input datasets to not overlap in
%       time for each unique DSI separately. Can therefore not handle
%       SOLO_L2_RPW-LFR-SBM1-CWF-E.
%       NOTE: Implemented via
%             bicas.tools.batch.autocreate_input_BPCIs()
%             calling
%             bicas.tools.batch.autocreate_many_BPCIs().
%   PROPOSAL: Define (hardcoded) list of DSIs "main" input datasets.
%       There is one main DSI in every SWM (one main input dataset in evey
%       BPCI). Separately for every SWM, iterate over all main datasets. For
%       every main dataset, find those other datasets (with input DSIs in
%       the SWM) which overlap with the main dataset in time.
%       --
%       PROPOSAL: Same list as INPUT_DSI_FOR_OUTPUT_TIME in functions
%           bicas.tools.batch.get_BPCI_output_path2()
%           bicas.tools.batch.default_get_BPCI_output_filename()
%           (latter to be phased out).
%   PROPOSAL: Given a set of input DSIs specified by a SWM, find any group
%             of datasets with some shared overlap.
%       NOTE: Datasets may be used in multiple BPCIs (for same SWM). This is
%             legitimately expected for
%             (1) CURRENT datasets,
%             (2) non-SBM datasets in SBM SWMs.
%           NOTE: Since L1R SBM1s overlap in time, there should legitimately
%                 be L2 SBM1s overlapping in time.
%       PRO: Symmetric w.r.t. datasets. Does not need to define any "main"
%            dataset for every SWM.
%       CON: Having two DSIs with time overlapping datasets in the same SWM
%            leads to processing four different combinations of datasets for
%            producing data for the same time interval.
%           Ex: Extreme but clear case:
%               Datasets 1a and 1b: Almost same time interval. Same DSI_1.
%               Datasets 2a and 2b: Almost same time interval. Same DSI_2.
%               ==> Output dataset(s)
%                 1a+2a ==> 3a
%                 1a+2b ==> 3b
%                 1b+2a ==> 3c
%                 1b+2b ==> 3d
%                 If DSI_1 determines the output filename, then
%                 3a and 3b will have the same filename (except version),
%                 3c and 3d will have the same filename (except version).
%           CON: This should never happen with real datasets. SBM1 is the
%                only dataset known to be overlapping with itself. CURRENT
%                datasets overlap with many other datasets because they are
%                long, not because they overlap with themselves.
%       CON: Unclear what the algorithm should be.
%
% PROPOSAL: Outsource implementation to new function "main2". "main" should only
%           pass on all arguments to "main2" and set
%           Bpa=bicas.tools.batch.BicasProcessingAccessImpl().
%   PRO: Good for automated tests. Tests for
%        bicas.tools.batch.run_BICAS_all_passes() could be converted into tests
%        for "main2".


% IMPLEMENTATION NOTE: bicas.tools.batch.main() calls
% bicas.create_default_BSO() directly which in turn uses irfu-matlab
% code (+irf/). This happens before BICAS is called, which itself
% initializes irfu-matlab paths using the same command, but then it is too
% late.
irf('check_path')

%==========
% Settings
%==========
DEFAULT_SETTINGS = [];
% How many days to extend the time coverage of CURRENT datasets after last
% timestamp (bias current setting).
DEFAULT_SETTINGS.currentDatasetExtensionDays = 0;
%------------------------------------------------------------------------
% BICAS settings
% --------------
% NOTE: These settings are always set in this code, thus overwriting any
% setting in config file.
% NOTE: Capitalized since that is the BICAS naming convention for
% BICAS-internal settings.
%------------------------------------------------------------------------
%DEFAULT_SETTINGS.bicasSetting_SWM_L1_L2_ENABLED   = 1;
%DEFAULT_SETTINGS.bicasSetting_SWM_L2_L3_ENABLED   = 1;
%
Settings = irf.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
clear DEFAULT_SETTINGS
assert(isnumeric(Settings.currentDatasetExtensionDays))



%======================================
% ASSERTIONS: Arguments (non-settings)
%======================================
% Useful to check existence of config file first since a faulty path will
% oterwise be found first much later.
irf.assert.file_exists(bicasConfigFile)
assert(isscalar(outputIsCdag) & islogical(outputIsCdag))
irf.assert.dir_exists(outputDir)
assert(iscell(inputPathsCa))



SwmArray = get_SWMs(bicasConfigFile);

bicasSettingsArgsCa = {};
%     bicasSettingsArgsCa(end+1:end+3) = {'--set', 'SWM.L1-L2_ENABLED',         sprintf('%i', Settings.bicasSetting_SWM_L1_L2_ENABLED)};
%     bicasSettingsArgsCa(end+1:end+3) = {'--set', 'SWM.L2-L3_ENABLED',         sprintf('%i', Settings.bicasSetting_SWM_L2_L3_ENABLED)};



switch(modeStr)
  case 'V01'
    fnVerAlgorithm = 'HIGHEST_USED';

    % IMPLEMENTATION NOTE: The filename version algorithm is more
    % generic. This code only creates datasets for which there is no
    % pre-existing dataset (ref.dir.). ==> All created datasets are V01,
    % when using this filename version algorithm. ==> Name "V01".

  case 'NEW_VERSION'
    fnVerAlgorithm = 'ABOVE_HIGHEST_USED';

  otherwise
    error('Illegal argument modeStr="%s".', modeStr)
end



% 31 = Format: 2021-03-19 20:08:34
fprintf('%s: Starting BICAS passes.\n', datestr(now, 31))
t = tic();

%====================
% CALL MAIN FUNCTION
%====================
Bpa = bicas.tools.batch.BicasProcessingAccessImpl();

BpcsArray = bicas.tools.batch.run_BICAS_all_passes(...
  Bpa, bicasSettingsArgsCa, ...
  bicasConfigFile, outputDir, referenceDir, inputPathsCa, ...
  fnVerAlgorithm, outputIsCdag, SwmArray, Settings);

%=======================
% Log BICAS error codes
%=======================
becArray = [BpcsArray.errorCode];
nNonError = nnz(becArray == 0);
nError    = nnz(becArray ~= 0);
fprintf('#BICAS calls, no error: %i\n', nNonError);
fprintf('              error:    %i\n', nError);

%================
% Log time usage
%================
nTpd = sum(arrayfun(@(Bpcs) numel(Bpcs.Bpci.outputsArray), BpcsArray));
wallTimeSec = toc(t);
fprintf('%s: Finished BICAS passes.\n',            datestr(now, 31));
fprintf('SPEED: Wall time: %.3f [h] = %.0f [s]\n', wallTimeSec / 3600, wallTimeSec);
fprintf('SPEED:            %.2f [s/TPD]\n',        wallTimeSec / nTpd);



if nError > 0
  error('%i of %i BICAS calls returned error.', ...
    nError, numel(becArray))
end
end



% Obtain information about BICAS s/w modes (SWMs) from BICAS using unofficial
% BICAS "interface" (in reality: Use knowledge of BICAS's implementation).
%
function SwmArray = get_SWMs(bicasConfigFile)
BSO = bicas.create_default_BSO();

% ID string used to inform BICAS SETTINGS of who set the setting. Only
% relevant for inspecting logs.
% NOTE: Exact string not really important.
% bicasSettingsSource = mfilename('fullpath');
%
%     BSO.override_value('SWM.L1-L2_ENABLED', ...
%         Settings.bicasSetting_SWM_L1_L2_ENABLED, ...
%         bicasSettingsSource)
%     BSO.override_value('SWM.L2-L3_ENABLED', ...
%         Settings.bicasSetting_SWM_L2_L3_ENABLED, ...
%         bicasSettingsSource)

bicas.override_settings_from_config_file(...
  bicasConfigFile, BSO, bicas.Logger('none', false));

BSO.make_read_only();

% NOTE: Converting SWML to array of SWMs.
SwmArray = bicas.swm.get_SWML(BSO).List;
end
