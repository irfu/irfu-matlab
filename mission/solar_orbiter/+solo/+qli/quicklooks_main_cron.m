%
% Wrapper around solo.qli.quicklooks_main() intended for being used by being
% called from system scripts (e.g. bash) for the purpose of cron jobs on
% brain/spis. The arguments have also been designed for this purpose and
% therefore all strings.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.quicklooks_main() instead.
%
%
% ARGUMENTS
% =========
% logoPath
%       Path to IRF logo image.
%       Normally located in irfu-matlab:
%       irfu-matlab/mission/solar_orbiter/+solo/irf_logo.png
%       Empty ==> Do not plot any logo.
% vhtDataDir
%       Path to directory containing VHT (velocity) .mat files.
% outputDir
%       Plots will be placed in subdirectories under this directory.
%       NOTE: Will create subdirectories if not pre-existing.
% runNonweeklyPlots, runWeeklyPlots
%       NOTE: STRINGS.
%       Whether to run ("1") or not run ("0") the resp. groups of plots.
% utcBegin, utcEnd : Strings.
%       Defines time interval for which quicklooks should be generated.
%
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
% * Happened when nobody worked (Sunday morning), i.e. low access hours.
% * Happened only 1h33m within each other's crashes.
% * Can be explained by disk error.
% * EXAMPLE 1 & 2: Could be re-run without triggering error for the same data.
% There are earlier reasons to believe that the nas24 disks are not always
% accessible (or not accessible quickly enough?) during low access hours,
% presumably due to being unmounted due to automounting.
% /Erik P G Johansson 2022-09-20
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function quicklooks_main_cron(...
  logoPath, vhtDataDir, outputDir, ...
  runNonweeklyPlots, runWeeklyPlots, utcBegin, utcEnd)

runNonweeklyPlots = interpret_argument_flag(runNonweeklyPlots);
runWeeklyPlots    = interpret_argument_flag(runWeeklyPlots);

% IMPLEMENTATION NOTE: Needed to make "DB" work. Necessary when calling from
% bash.
irf

%===============================================================
% Configure Solar Orbiter database from which data will be used
%===============================================================
% NOTE: System-dependent configuration!
solo.db_init('local_file_db', '/data/solo/');
solo.db_init('local_file_db', '/data/solo/data_irfu');
% Setup cache
solo.db_init('db_cache_size_max', 4096)
solo.db_cache('on', 'save')

%======
% Plot
%======
solo.qli.quicklooks_main(...
  logoPath, vhtDataDir, outputDir, ...
  runNonweeklyPlots, runWeeklyPlots, utcBegin, utcEnd)
end



% Interpret argument for main function interface. Intended accept and normalize
% arguments which are either
% (1) MATLAB-friendly (numeric/logical), or
% (2) bash script-friendly (strings).
%
function value = interpret_argument_flag(arg)
assert(isscalar(arg), 'Flag argument is not scalar.')

if ischar(arg) && arg=='0'
  value = false;
elseif ischar(arg) && arg=='1'
  value = true;
else
  error('Can not interpret argument flag. Illegal format.')
end
end
