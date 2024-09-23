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
%   MATLAB). All arguments are therefore strings.
% * See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT for the meaning
%   of abbreviations.
%
%
% CONFIG FILE
% ===========
% One of the function arguments refers to a JSON config file which stores
% system-dependent values. Below example contains comments which have to be
% removed before using it as a file (comments are illegal in JSON files(!)):
% ------------------------------------------------------------------------------
% {
%   # Path to IRF logo or null (not quoted string) if no logo.
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
%   # It should be enough to point to the 24h quicklooks subdirectory. Directory
%   # is searched recursively and uses filenames to recognize files.
%   # This parameter does not specify the output directory, though the value may
%   # refer to the same directory as the output directory.
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
% ------------------------------------------------------------------------------
%
%
% ARGUMENTS
% =========
% NOTE: All arguments are strings.
% --
% outputDir
%       Path to output directory.
% configFilePath
%       Path to config file (see other section in documentation).
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       One-character strings. Whether to generate ("1" or "0") non-weekly (2h,
%       6h, 24h) quicklooks and/or weekly quicklooks.
%       NOTE: STRINGS (more natural when calling from bash script).
% operationId
%       String constant which specifies what to do with list of dates, once
%       obtained. Either
%       "LIST":
%           Only list the dates for which quicklooks would be generated if other
%           arguments are identical.
%       "GENERATE":
%           Generate actual quicklooks for dates indirectly specified by
%           arguments "dasaid" and "varargin".
% dasaid
%       String constant which specifies which DASA (algorithm) should be used
%       for obtaining dates.
% varargin
%       Arguments used by the selected DASA (dasaid). Meaning varies depending
%       on "dasaid". See section below.
%
%
% SYNTAX FOR ARGUMENTS "dasaid" AND varargin
% ==========================================
% The combination of arguments "dasaid" and varargin specifies how to generate
% the list of dates for which quicklooks should be generated. "dasaid" specifies
% the "algorithm" (DASA) and varargin specifies the arguments for that specific
% algorithm. All arguments are strings. Dates are on the format YYYY-MM-DD. When
% specifying a range of dates, then the first date is inclusive and the last
% date is exclusive.
% --
% dasaid              varargin
% -----------------------------------------------------------------------------
% "TIME_INTERVAL"     firstDate  lastDate
%     Continuous time interval.
% "LOGS"              logType1 ... logTypeN
%     Search log file(s) for dataset filenames.
%     LOG_TYPE_i = String constant specifying type of log. "LESIA" or "SOAR".
%     Corresponds to keys into solo.qli.batch.const.SOURCE_DSI_DICT.
% "DMRQ"              maxNbrOfDays  firstDate  lastDate
%     Use dates for which the file modification dates of the newest dataset is
%     more recent than the file modification date of the oldest corresponding
%     quicklook.
% "QLI_FMD_INTERVAL"  maxNbrOfDays  oldestFmd  newestFmd
%     Use dates within a specific interval(!) of file modification dates (FMD).
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
if ~isempty(Config.irfLogoPath)
  irf.assert.file_exists(Config.irfLogoPath)
end



%=============================================================
% Configure SolO DB from which data will be used for plotting
%=============================================================
% NOTE: brain/spis: solo.db_init('local_file_db', '/data/solo/'); implies only
% using /data/solo/remote/data/ and /data/solo/soar/.
solo.db_init('local_file_db', Config.soloDbDirPath);

% Setup cache
solo.db_init('db_cache_size_max', 4096)   % Unit: MiB.
solo.db_cache('on', 'save')



%=============================
% Configure the parallel pool
%=============================
% NOTE: This requires there to not be a pre-existing parallel pool (will
% otherwise crash).
% --

% Delete any optionally pre-existing parallel pool.
delete(gcp('nocreate'))

% spmdEnabled=false : EXPERIMENTAL. May prevent solo.qli.batch from crashing on
% anna.irfu.se. /2024-08-16.
% NumWorkers=2 : EXPERIMENTAL. May prevent solo.qli.batch from crashing on
% anna.irfu.se. /2024-08-19.
% NOTE: spmdEnabled=false makes it illegal to use the "spmd" commands.
parpool('SpmdEnabled', false);
% if isunix()
%   [~, hostName] = system('hostname');
%   hostName      = strtrim(hostName);
%   isServerAnna  = strcmp(hostName, 'anna');
% else
%   isServerAnna  = false;
% end
% ppArgsCa = {'Processes', 'SpmdEnabled', false};
% if isServerAnna
%   % Set number of workers to a lower value than the default on anna (4).
%   ppArgsCa = {'Processes', 1, 'SpmdEnabled', false};
% else
%   ppArgsCa = {'Processes', 'SpmdEnabled', false};
% end
% parpool(ppArgsCa{:});



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
