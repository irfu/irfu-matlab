%
% Wrapper around solo.qli.generate_quicklooks_batch() intended for being run
% on brain/spis for the purpose of cron jobs and official generation of
% quicklooks.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT.
%
% NOTE: This function is NOT intended to be called from bash/the OS. There are
%       other MATLAB scripts for that.
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory.
%       IMPLEMENTATION NOTE: Not hardcoding this is useful for manual testing of
%       the setup and for manual large-scale processing.
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       Scalar logical. Whether to generate non-weekly (2h, 6h, 24h) quicklooks
%       and/or weekly quicklooks.
% DaysDtArray
%       datetime column array. Array of UTC midnights representing the beginning
%       of days for which to generate quicklooks.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function generate_quicklooks(...
  outputDir, generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
%
% NOTES ON CRASHES ON SPIS & BRAIN
% ================================
% Long runs on brain/spis may crash due to not being able to access disk
% (nas24), which is unrelated to the QLI code itself.
%
% EXAMPLE 1: Crash 2022-09-18 (Sunday) 06:57:50, brain:
% Crashed while processing data after 2022-03-23.
% """"""""
% purging solo_L3_rpw-bia-efield-10-seconds_20220313_V01.cdf
% purging solo_L3_rpw-bia-density-10-seconds_20220313_V01.cdf
% purging solo_L2_swa-pas-grnd-mom_20220313_V02.cdf
% purging solo_L2_rpw-tnr-surv-cdag_20220313_V05.cdf
% exception.message=No cdf files specified
% dataobj, row 72
% rcdf, row 229
% read_TNR, row 53
% quicklooks_24_6_2_h, row 203
% quicklooks_main, row 209
% quicklooks_main_cron, row 54
% Command exited with non-zero status 99
% """"""""
%
% EXAMPLE 2: Crash 2022-09-18 (Sunday) 08:30:13, spis
% Crashed while processing data after 2021-12-17.
% """"""""
% purging solo_L3_rpw-bia-efield-10-seconds_20211209_V01.cdf
% purging solo_L3_rpw-bia-density-10-seconds_20211209_V01.cdf
% purging solo_L2_swa-pas-grnd-mom_20211208_V03.cdf
% purging solo_L2_swa-pas-eflux_20211208_V02.cdf
% purging solo_L2_rpw-tnr-surv-cdag_20211208_V04.cdf
%   [warning: solo_db.get_variable/append_sci_var(131)] Discarded 3745 data points
% exception.message=SPICE(INVALIDVALUE): [spkpos_c->SPKPOS->SPKEZP->SPKAPO->SPKGPS->SPKPVN->SPKR19] Window size in type 19 segment was 0; must be in the range 2:14 for subtype 0. Mini-segment index is 113. Failure occurred at input vector index 7. (CSPICE_N0067)
% cspice_spkpos, row 630
% get_position, row 93
% get_SolO_pos, row 320
% quicklooks_main, row 189
% quicklooks_main_cron, row 54
% Command exited with non-zero status 99
% """"""""
%
% NOTE: Above crashes...
% * Happened after 300h+ (12-13 days) of execution.
% * Happened when nobody else was working (Sunday morning), i.e. low-access
%   hours.
% * Happened only 1h33m within each other's crashes.
% * Can be explained by disk error.
% * EXAMPLE 1 & 2: Could be re-run without triggering error for the same data.
% There are earlier reasons to believe that the nas24 disks are not always
% accessible (or not accessible quickly enough?) during low access hours,
% presumably due to being unmounted due to automounting.
% /Erik P G Johansson 2022-09-20


VHT_DIR = '/data/solo/data_yuri/';

% Path to IRF logo, relative to the irfu-matlab root.
% NOTE: The IRF logo is not part of the irfu-matlab git repo, but this code still
% requires it to be located inside the corresponding directory.
IRF_LOGO_RPATH = 'mission/solar_orbiter/+solo/+qli/+offgen/irf_logo.png';



assert(islogical(generateNonweeklyQuicklooks))
assert(islogical(generateWeeklyQuicklooks))



irfLogoPath = fullfile(irf('path'), IRF_LOGO_RPATH);
irf.assert.file_exists(irfLogoPath)

% IMPLEMENTATION NOTE: Needed to make "DB" work. Necessary when calling from
% bash.
irf()

% Set log level. Must use "notice" to not miss relevant log messages, in
% particular which time interval is currently being plotted (which is useful for
% debugging crashes during long batch runs).
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
solo.qli.generate_quicklooks_batch(...
  irfLogoPath, VHT_DIR, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
end
