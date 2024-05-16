%
% Wrapper around solo.qli.batch.generate_quicklooks_syntax() intended for
% being run on brain/spis for the purpose of OFFICIAL GENERATION of quicklooks,
% e.g. in cron jobs.
%
%
% NOTES
% =====
% * This script is NOT intended to be called from MATLAB by the average
%   user as is. See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT.
%   It may however be seen as an example of how to call
%   solo.qli.batch.generate_quicklooks_syntax().
% * This script is intended for being used for official official generation of
%   QLIs, including cron jobs.
% * This function is designed for being called from bash/the OS (not from
%   MATLAB).
% * See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT for the meaning
%   of abbreviations.
%
%
% CONFIG FILE
% ===========
% JSON file which stores system-dependent values.
% Example, plus comments which are illegal in JSON files(!):
% {
%   # Path to IRF logo.
%   # See solo.qli.batch.generate_quicklooks().
%   "irfLogoPath": "/home/erjo/so_irfu-matlab_qli_IRFpl/mission/solar_orbiter/+solo/+qli/+batch/irf_logo.png",
%
%   # Path to VHT directory.
%   # See solo.qli.batch.generate_quicklooks().
%   "vhtDir":       "/data/solo/data_yuri/",
%
%   # Path used by "SolO DB" (submitted to solo.db_init('local_file_db', ...);
%   "solo_db_init_local_file_db": "/data/solo/",
%
%   # Directories with datasets which will effectively be used by
%   # solo.db_get_ts() and solo.db_list_files() for locating requested data.
%   # Must be consistent with "solo_db_init_local_file_db".
%   # See solo.qli.batch.interface.get_days_from_DMRQ().
%   "datasetDirs": [
%     "/data/solo/remote/data",
%     "/data/solo/soar"
%   ],
%
%   # Directory with quicklooks (QLIs) which is used for determining the
%   # file-modification dates (FMDs) of pre-existing QLIs.
%   # Should be enough to point to the 24h quicklooks.
%   # This is NOT the output directory!
%   "fmdQliDir":    "/data/juice/EJ_temp/quicklooks_SOLAR_ORBITER/www/24h",
%
%   # Paths with globbing filename patterns for log files which contain
%   # the names of new datasets.
%   # See solo.qli.batch.interface.get_days_from_logs().
%   "logFileDirPatterns": {
%     "LESIA": "/home/erjo/logs/pull.so.data.cron.brain.*.log",
%     "SOAR":  "/home/erjo/logs/so_soar_irfu_mirror_sync.*.log"
%   }
% }
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory.
% configFilePath
%       Path to config file (see other section in documentaiton).
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       One-character strings. Whether to generate ("1" or "0") non-weekly (2h,
%       6h, 24h) quicklooks and/or weekly quicklooks.
%       NOTE: STRINGS.
% operationId
%       String constant which specifies what to do with list of dates, once
%       obtained. Either "LIST" or "GENERATE".
% dasaid
%       String constant which specifies which DASA (algorithm) should be used
%       for obtaining dates.
% varargin
%       Arguments used by the selected DASA (dasaid). Meaning varies depending
%       on "dasaid".
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_shell(...
  outputDir, configFilePath, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  operationId, dasaid, varargin)

% IMPLEMENTATION NOTE: Argument outputDir comes before configFilePath since this
% makes it easier to call MATLAB code from so_qli2.

% IMPLEMENTATION NOTE: Needed to make "SolO DB" work. Necessary when calling
% from bash.
irf()

% Set log level. Must use "notice" to not miss relevant non-error/warning log
% messages, in particular which time interval is currently being plotted which
% is useful for debugging crashes during long batch runs.
irf.log('notice')



generateNonweeklyQuicklooks = solo.qli.batch.interface.interpret_boolean_flag(generateNonweeklyQuicklooks);
generateWeeklyQuicklooks    = solo.qli.batch.interface.interpret_boolean_flag(generateWeeklyQuicklooks);
dasaArgumentsCa             = varargin(:);    % Column array.

Config = solo.qli.batch.utils.read_config_file(configFilePath);
irf.assert.file_exists(Config.irfLogoPath)



%=============================================================
% Configure SolO DB from which data will be used for plotting
%=============================================================
% NOTE: brain/spis: solo.db_init('local_file_db', '/data/solo/'); implies only
% using /data/solo/remote/data/ and /data/solo/soar/.
solo.db_init('local_file_db', Config.soloDbDirPath);

% Setup cache
solo.db_init('db_cache_size_max', 4096)   % Unit: MiB.
solo.db_cache('on', 'save')



%=====================
% Generate quicklooks
%=====================
solo.qli.batch.generate_quicklooks_syntax(...
  Config, ...
  solo.qli.batch.GenerateQuicklooksImplementation(), ...
  solo.qli.batch.FileSystemReaderImplementation(), ...
  outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  operationId, dasaid, dasaArgumentsCa)
end
