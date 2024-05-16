%
% Generate multiple types of quicklooks (24 h, 6 h, 2 h, weekly; files) for an
% array of explicitly specified dates.
%
%
% NOTES
% =====
% * Requires solo.db_init() to have been properly used to initialize dataset
%   lookup before calling this function.
% * Uses SPICE implicitly, and therefore relies on some path convention used by
%   irfu-matlab for where to find SPICE kernels. Not sure which convention, but
%   presumably it does at least find /data/solo/SPICE/.
% * Creates subdirectories to the output directory if not pre-existing.
% * Note: "Weekly"/7-day quicklooks always begin with a specific hardcoded
%   weekday (Wednesday as of 2024-03-22).
% * Overwrites pre-existing quicklook files without warning.
% * ~BUG: SolO DB (solo.db_get_ts() etc.) requires the caller to specify
%   dataset_ID plus "-cdag" if present, but does not raise exception if wrong.
%   ==> The code requires the datasets searched by SolO DB to have/lack "-cdag"
%   exactly as specified (hardcoded). If they are not, then data appears to not
%   be present and no exception is raised! This means that the code might not
%   recognize datasets for the path specified with solo.db_init().
% * The code uses irf.log(). The caller may want to set the log level before
%   calling the code.
%
%
% ARGUMENTS
% =========
% irfLogoPath
%       Path to IRF logo. Empty ==> Do not plot any logo.
%       IMPLEMENTATION NOTE: The IRF logo should only be used for official
%       quicklooks.
% vhtDataDir
%       Path to directory containing VHT (velocity) data files
%       V_RPW_1h.mat and V_RPW.mat. Typically brain:/data/solo/data_yuri/
% outputDir
%       Quicklooks will be placed in subdirectories under this directory.
%       NOTE: Will create subdirectories if not pre-existing.
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       Scalar logical. Whether to generate non-weekly (2h, 6h, 24h) quicklooks
%       and/or weekly quicklooks.
%       Useful for testing and not re-running unnecessary time-consuming
%       quicklooks.
% DaysDtArray
%       Column array of datetime objects. UTC, midnight only
%       (TimeZone = "UTCLeapSeconds").
%       Defines days for which quicklooks should be generated. A day is
%       specified by the midnight which begins that day.
%       Weekly quicklooks (if enabled) will be produced for all weeks which
%       contain at least one day specified here.
% Gql
%       solo.qli.batch.GenerateQuicklooks* object. Should be an instance of
%       solo.qli.batch.GenerateQuicklooksImplementation object in the nominal
%       case. Argument exists to facilitate automated tests.
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function generate_quicklooks(...
  irfLogoPath, vhtDataDir, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray, Gql)

% NOTE: Data begins on 2020-02-12=Wednesday.
% ==> There is no SPICE data on Monday-Tuesday before this date.
% ==> Code fails for week Monday-to-Sunday.
%     PROPOSAL: Additionally round start time up to start date.
%
% PROPOSAL: Some way of handling disk access error?
%   PROPOSAL: try-catch plot code once (weekly or non-weekly quicklooks function).
%             Then try without catch a second time, maybe after delay.
%             If the first call fails due to disk access error, it might still
%             trigger automount which makes the second attempt succeed.
%   PROPOSAL: Before loading data with SolO DB, use same trick as in bash
%             scripts for triggering automounting of NAS before every time
%             SolO DB is called. -- IMPLEMENTED BUT UNTESTED /2024-03-07
%     NOTE: Has seen presumed automount fail errors in long quicklook generation
%           runs.
%         Ex: so_qli.2024-02-27_20.11.02.log
%
% PROPOSAL: Parallelize loop over days and weeks.
%
% PROPOSAL: Always call solo.db_get_ts() both with and without "-cdag", to make
%           sure that the code does select non-existing datasets.
%           Use pre-existing solo.qli.utils.db_get_ts() wrapper.
%
% PROPOSAL: Replace argument vhtDataDir --> vhtData1hFilePath,
%           vhtData6hFilePath.
%   PRO: Inconsistent to have hardcoded filenames but argument for directory.
%   PRO: cron setup should specify filenames.
%   PRO: generate_quicklooks_*_using_DB_SPICE() alraedy have arguments for paths
%        to he files directly.
%
% PROPOSAL: Log data date of re-thrown exception (previously caught exception).
%
% PROBLEM: How handle that simultaneous batch processing runs (using potentially
%          different methods for deriving list of dates)?
%   PROPOSAL: Create temporary file in standard location which informs other
%             processes of dates for which the current process will generate quicklooks.
%             When a processes begins, it reads pre-existing temporary file(s)
%             and removes those from its own list of dates to generate.
%   PROPOSAL: Check FMD for quicklooks before generating them. If FMD is more
%             recent than the time current process launched, then do not
%             generate for that date.
%     CON: Only works easily for 24h quicklooks (which are not always generated).
%   PROPOSAL: Have some kind of queue system: One queue of dates.
%
% PROPOSAL: Before starting generation of quicklook(s) for each day/week, check
%           timestamp of output quicklook(s). Only actually generate quicklook
%           if the pre-existing quicklooks have an FMD older than a specified
%           reference timestamp (probably) equal to the beginning of the function call
%           (solo.qli.batch.generate_quicklooks).
%   PRO: Simultaneous batch runs will not generate the same quicklooks (mostly).
%     CON: Not true. The scheme will alleviate but not solve problem of
%          simultaneously running batch runs.
%       Ex: Two batch runs A & B generate for the same dates. A starts before B.
%           A will produce some quicklooks before B launches. These will be
%           regenerated by B. When A & B runs simultaneously, they will not
%           regenerate quicklooks created by the other run.
%   PROPOSAL: Set reference timestamp to FMD of relevant log(s).
%   PROPOSAL: Caller sets the reference timestamp.
%
%
% generate_quicklooks_24h_6h_2h(), generate_quicklook_7day()
% ============================================================
% PROBLEM: Lots of cases of, and checks for data that may or may not be missing.
% PROBLEM: Missing data can be represented as [] (0x0), rather than e.g. Nx0.
%   Stems from solo.db_get_ts() (?) not finding any data for selected time
%   interval. Can not know the dimensions of missing data if can not find any
%   CDF at all.
%
%
% Examples of missing data
% ========================
% POLICY: Incomplete list of examples of missing data that could be useful for
%         testing(?).
% POLICY: Not necessarily the entire data gaps. Only segments that trigger
%         special handling or errors.
% NOTE: SolO data gaps may be filled laters.
% --
% 2022-03-12T07:22:03.090373000Z -- 2022-03-13T00
%   Missing data.Tpas



%=========================
% Assertions on arguments
%=========================
assert(...
  islogical(generateNonweeklyQuicklooks) & isscalar(generateNonweeklyQuicklooks), ...
  'Argument generateNonweeklyQuicklooks is not a scalar logical.')
assert(...
  islogical(generateWeeklyQuicklooks   ) & isscalar(generateWeeklyQuicklooks), ...
  'Argument generateWeeklyQuicklooks is not a scalar logical.')
assert(iscolumn(DaysDtArray), 'Argument DaysDtArray is not column array (Nx1).')
solo.qli.utils.assert_UMD_DT(DaysDtArray)
assert(isa(Gql, 'solo.qli.batch.GenerateQuicklooksAbstract'))



% "Normalize"
DaysDtArray = unique(DaysDtArray);



%==============================
% Miscellaneous initialization
%==============================
OutputPaths = solo.qli.utils.create_output_directories(outputDir);

automountTriggeringPathsCa = {vhtDataDir, outputDir};
if ~isempty(irfLogoPath)
  automountTriggeringPathsCa{end+1} = irfLogoPath;
end



% Log arguments
irf.log('n', sprintf('irfLogoPath                   = "%s"', irfLogoPath))
irf.log('n', sprintf('vhtDataDir                    = "%s"', vhtDataDir))
irf.log('n', sprintf('outputDir                     = "%s"', outputDir))
irf.log('n', sprintf('generateNonweeklyQuicklooks   = %d',   generateNonweeklyQuicklooks))
irf.log('n', sprintf('generateWeeklyQuicklooks      = %d',   generateWeeklyQuicklooks))
irf.log('n', sprintf('numel(DaysDtArray)            = %d',   numel(DaysDtArray)))
% Log misc. variables
irf.log('n', sprintf('automountTriggeringPathsCa    = {%s}', strjoin(automountTriggeringPathsCa, ', ')))
% Log selected constants.
irf.log('n', sprintf('B_SPECTRA_ENABLED             = %d',   solo.qli.const.B_SPECTRA_ENABLED))
irf.log('n', sprintf('NONWEEKLY_ALL_PLOTS_ENABLED   = %d',   solo.qli.const.NONWEEKLY_ALL_PLOTS_ENABLED))
irf.log('n', sprintf('NONWEEKLY_6H_2H_PLOTS_ENABLED = %d',   solo.qli.const.NONWEEKLY_6H_2H_PLOTS_ENABLED))
% Log current working directory so that relative paths can be interpreted.
irf.log('n', sprintf('Current working directory     = %s',   pwd))



tSec = tic();

% Array of plotting exceptions caught.
PlotExcArray = MException.empty(1, 0);



%=============================================
% Run the code for 2-, 6-, 24-hour quicklooks
%=============================================
if generateNonweeklyQuicklooks
  % Daily time-intervals

  % This is the .mat file containing RPW speeds at 1h resolution.
  % The file should be in the current path. This file can be found in
  % brain:/solo/data/data_yuri/.
  vhtFile1hPath = fullfile(vhtDataDir, solo.qli.const.VHT_1H_DATA_FILENAME);

  for iDay = 1:length(DaysDtArray)
    DayDt = DaysDtArray(iDay);

    try
      trigger_automounting(automountTriggeringPathsCa)

      irf.log('n',         '============================================================')
      irf.log('n', sprintf('Calling 24h/6h/2h plot function for %s', string(DayDt)))
      irf.log('n',         '============================================================')
      tBeginSec = tic();
      Gql.generate_quicklooks_24h_6h_2h_using_DB_SPICE(DayDt, vhtFile1hPath, OutputPaths, irfLogoPath)
      solo.qli.utils.log_time('Time to generate one day''s 24h/6h/2h quicklooks', tBeginSec);
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(Exc)
    end
  end
end



%===================================
% Run the code for weekly overviews
%===================================
if generateWeeklyQuicklooks

  % Derive weeks from specified days (midnights which begin 7-day periods).
  WeeksDtArray = solo.qli.utils.derive_weeks(DaysDtArray, solo.qli.const.FIRST_DAY_OF_WEEK);

  % This is the .mat file containing RPW speeds at 6h resolution.
  % The file should be in the same folder as this script (quicklook_main).
  vhtFile6hPath = fullfile(vhtDataDir, solo.qli.const.VHT_6H_DATA_FILENAME);

  for iWeek = 1:numel(WeeksDtArray)
    WeekDt = WeeksDtArray(iWeek);

    try
      trigger_automounting(automountTriggeringPathsCa)

      irf.log('n',         '========================================================')
      irf.log('n', sprintf('Calling 7-day plot function for %s', string(WeekDt)))
      irf.log('n',         '========================================================')
      tBeginSec = tic();
      Gql.generate_quicklook_7days_using_DB_SPICE(WeekDt, vhtFile6hPath, OutputPaths.dir1w, irfLogoPath)
      solo.qli.utils.log_time('Time to generate one 7-day quicklook', tBeginSec);
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(Exc)
    end
  end
end



wallTimeSec     = toc(tSec);
wallTimeHours   = wallTimeSec/3600;
nQuicklooksDays = numel(DaysDtArray);

% NOTE: Execution speed may vary by orders of magnitude depending on settings
% (nonweekly vs weekly plots). May therefore want scientific notation.
irf.log('n', sprintf('Wall time used:                       %g [h] = %g [s]', wallTimeHours, wallTimeSec));
irf.log('n', sprintf('Wall time used per day of quicklooks: %g [h/day]',      wallTimeHours / nQuicklooksDays));



if solo.qli.const.CATCH_PLOT_EXCEPTIONS_ENABLED && ~isempty(PlotExcArray)
  %fprintf(2, 'Caught %i plotting exceptions.\n', numel(PlotExcArray))
  %fprintf(2, 'Rethrowing old (last) exception.\n')
  irf.log('c', sprintf('Caught %i plotting exceptions.', numel(PlotExcArray)))
  irf.log('c', 'Rethrowing old (last) exception.')

  % NOTE: This does display (stderr) the stack trace for position
  % of the *ORIGINAL* error.
  rethrow(PlotExcArray(end))
end

end    % function



% Handle *PLOTTING* exception.
%
% Historically, the plotting code has caused many exceptions. One may want
% different behaviour depending on context.
%
% Production: Produce as many plots as possible.
%             => Catch exception and continue.
% Testing:    Crash on first exception so that it can be fixed.
%             => Rethrow exception as soon as possible.
function handle_plot_exception(Exc)
if solo.qli.const.CATCH_PLOT_EXCEPTIONS_ENABLED
  % Print stack trace without rethrowing exception.
  % One wants that in the log.
  % NOTE: fprintf(FID=2) => stderr
  %fprintf(2, 'Caught plotting error without rethrowing it.\n')
  %fprintf(2, 'Plot error/exception: "%s"\n', Exc.message)
  irf.log('w', 'Caught plotting error without rethrowing it.')
  irf.log('w', sprintf('Plot error/exception: "%s"', Exc.message))

  for i = 1:numel(Exc.stack)
    s = Exc.stack(i);
    %fprintf(2, '    Error in %s (line %i)\n', s.name, s.line)
    irf.log(...
      solo.qli.const.LOG_LEVEL_CAUGHT_EXCEPTIONS, ...
      sprintf('    Error in %s (line %i)', s.name, s.line))
  end
else
  rethrow(Exc)
end
end



% Try to trigger automounting of disks in order to avoid file-reading problems
% for long batch runs.
%
%
% CONTEXT
% =======
% The automounting which should make IRFU's /data/solo/ available on brain/spis
% (IRFU) does not always work as intended. It appears to not always work fast
% enough for reading files on rare occassions, leading to apparent
% can-not-read-file errors for existing files when making long runs, in
% particular on out-of-office hours(?). Other code which has faced this problem
% has seemingly resolved it by using Linux commands to list the directory
% contents and ignore the result before using /data/solo/ (such as "ls
% /data/solo/ >> /dev/null").
%
function trigger_automounting(pathsCa)
for i = 1:numel(pathsCa)
  path = pathsCa{i};

  irf.log('n', sprintf(...
    'Trying to trigger automounting, if not already mounted: %s', path ...
    ))
  junk = dir(path);
end
end