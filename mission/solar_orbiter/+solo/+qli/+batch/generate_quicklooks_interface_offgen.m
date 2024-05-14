%
% Wrapper around solo.qli.batch.generate_quicklooks_interface() intended for
% being run on brain/spis for the purpose of OFFICIAL GENERATION of quicklooks,
% e.g. in cron jobs.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user as is. See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT.
%       It may however be seen as an example of how to call
%       solo.qli.batch.generate_quicklooks_interface().
%
% NOTE: This script relies on hardcoded information about the system setup used
%       for official generation. This function is the intended store of
%       hardcoded information that is relevant ONLY to the OFFICIAL GENERATION.
%
% NOTE: This function is intended to be called from bash/the OS.
%
%
% ARGUMENTS
% =========
% varargin
%       Arguments which are passed on to
%       solo.qli.batch.generate_quicklooks_interface() (all but "Settings").
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_interface_offgen(varargin)
% PROPOSAL: Implement support for reading a config file (e.g. JSON) containing
%           the values hardcoded in this file. Then abolish the hardcoded values
%           (in the code).
%   PROPOSAL: Implement config file support in
%             solo.qli.batch.generate_quicklooks_interface()
%             and abolish this wrapper.
%     CON: Makes it harder to call
%           solo.qli.batch.generate_quicklooks_interface() from other MATLAB code.
%       CON: Is not meant to be called by other MATLAB code, except test code.
%       PRO: Settings.Gql object can not be stored in config file.
%         CON-PROPOSAL: Create wrapper which only sets GQL.
%           CON: Yet another file and level in the function call hierarchy.
%                <=> Can not abolish this function.

% IMPLEMENTATION NOTE: Not hardcoding outputPath since this is useful for
% manual testing of the setup and for manual large-scale processing.

% Path to IRF logo, relative to the irfu-matlab root.
% NOTE: The IRF logo is not part of the irfu-matlab git repo, but this code
% still requires it to be located inside the corresponding directory.
IRF_LOGO_RPATH  = 'mission/solar_orbiter/+solo/+qli/+batch/irf_logo.png';

VHT_DIR         = '/data/solo/data_yuri/';

% Nominal directory.
%FMD_QLI_DIR     = '/data/solo/quicklooks/www/24h';
% Directory used while the NAS8 SolO partition still has problems.
FMD_QLI_DIR     = '/data/juice/EJ_temp/quicklooks_SOLAR_ORBITER/www/24h';

% Directories with datasets which will be used by solo.db_get_ts() and
% solo.db_list_files().
DATASET_DIRS_CA = {
  '/data/solo/remote/data';    % LESIA
  '/data/solo/soar';           % SOAR
};

% Paths with globbing filename patterns for log files which contain the names of
% new datasets.
LOG_FILE_DIR_PATTERN_DICT          = dictionary();
LOG_FILE_DIR_PATTERN_DICT('LESIA') = '/home/erjo/logs/pull.so.data.cron.brain.*.log';
LOG_FILE_DIR_PATTERN_DICT('SOAR')  = '/home/erjo/logs/so_soar_irfu_mirror_sync.*.log';



Settings = [];
Settings.Gql                    = solo.qli.batch.GenerateQuicklooksImplementation();
Settings.Fsr                    = solo.qli.batch.FileSystemReaderImplementation();
Settings.vhtDir                 = VHT_DIR;
Settings.irfLogoPath            = fullfile(irf('path'), IRF_LOGO_RPATH);
Settings.LogFileDirPatternDict  = LOG_FILE_DIR_PATTERN_DICT;
Settings.datasetDirsCa          = DATASET_DIRS_CA;
Settings.fmdQliDir              = FMD_QLI_DIR;

irf.assert.file_exists(Settings.irfLogoPath)



% IMPLEMENTATION NOTE: Needed to make "DB" work. Necessary when calling from
% bash.
irf()

% Set log level. Must use "notice" to not miss relevant non-error/warning log
% messages, in particular which time interval is currently being plotted which
% is useful for debugging crashes during long batch runs.
irf.log('notice')

%===============================================================
% Configure Solar Orbiter database from which data will be used
%===============================================================
% NOTE: System-dependent configuration!
% IMPLEMENTATION NOTE: Only uses /data/solo/remote/data/ and /data/solo/soar/
% (and not /data/solo/data_irfu/) since data_irfu/ (1) has less reliable data
% ("bleeding edge"), and (2) is (somewhat) ~frequently reprocessed in large
% sets which would increase automatic QLI processing a lot.

solo.db_init('local_file_db', '/data/solo/');

% Setup cache
solo.db_init('db_cache_size_max', 4096)   % Unit: MiB.
solo.db_cache('on', 'save')



%=====================
% Generate quicklooks
%=====================
solo.qli.batch.generate_quicklooks_interface(Settings, varargin{:})
end
