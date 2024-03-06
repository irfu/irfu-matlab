%
% Generate multiple IRFU-local types of quicklooks (files) for SolO data (not
% just RPW) for a specified time interval.
%
% Intended for batch processing, e.g. being called from bash script, e.g. cron
% job via a MATLAB wrapper script.
%
%
% NOTES
% =====
% * Requires solo.db_init() to have been properly used to initialize dataset
%   lookup.
% * Uses SPICE implicitly, and therefore relies on some path convention. Not
%   sure which, but presumably it does at least find /data/solo/SPICE/.
% * Uses solo.read_TNR() indirectly which in turns relies on a hardcoded
%   path to "/data/solo/remote/data/L2/thr/" and selected subdirectories.
% * Creates subdirectories to the output directory if not pre-existing.
% * Note: 7-day quicklooks always begin with a specific hardcoded weekday
%   (Wednesday as of 2023-07-24).
% * Overwrites pre-existing quicklook files without warning.
%
%
% ARGUMENTS
% =========
% vhtDataDir
%       Path to directory containing VHT (velocity) .mat files
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
%       Column array of datetime objects. UTC, midnight only.
%       Defines days for which quicklooks should be generated. A day is
%       specified by the midnight which begins that day.
%       Weekly quicklooks (if enabled) will be produced for all weeks which
%       contain at least one day specified here.
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function quicklooks_main(...
  vhtDataDir, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
%
% NOTE: Mission begins on 2020-02-12=Wednesday.
% ==> There is no SPICE data on Monday-Tuesday before this date.
% ==> Code fails for week Monday-to-Sunday.
%     PROPOSAL: Additionally round start time up to start date.
%
% PROPOSAL: Better name: quicklooks_main?
%   NEED: Consistent with various wrappers.
%       ~main, ~plot, ~generate, ~qli, quicklooks
%       generate_quicklooks
%   PROPOSAL:
%       qli.generate_quicklooks()
%           This file.
%           Input: List of arbitrary dates.
%       qli.cron.*():
%           Wrapper around generate_quicklooks() which configures SolO DB for brain/spis.
%           Input: List of arbitrary dates.
%       qli.cron.*_bash_manual()
%           Wrapper around *() for being called from bash for manual
%           generation of quicklooks.
%           Input: Two timestamp arguments specifying one continuous range of dates.
%       qli.cron.*_bash_logs()
%       qli.cron.*_bash_from_logs()
%       qli.cron.*_bash_log_updates()
%       qli.cron.*_bash_from_log_updates()
%           Input: None?
%           PROPOSAL: Log file patterns?
%               CON: Varying set of log files.
% PROPOSAL: Rename
%       qli.quicklooks_24_6_2_h()
%       qli.quicklooks_7days()     # NOTE: Generates ONE quicklook.
%     NOTE: Files are used by JB.
%     qli.generate_quicklooks_*
%
% PROPOSAL: Some way of handling disk access error?
%   PROPOSAL: try-catch plot code once (weekly or non-weekly quicklooks function).
%             Then try without catch a second time, maybe after delay.
%             If the first call fails due to disk access error, it might still
%             trigger automount which makes the second attempt succeed.
%   PROPOSAL: Before loading data with SolO DB, use same trick as in bash
%             scripts for triggering automounting of NAS.
%     NOTE: Has seen presumed automount fail errors in long quicklook generation
%           runs.
%         Ex: so_qli.2024-02-27_20.11.02.log
%
% PROPOSAL: Log time consumption for each call to plot functions.
%   PROPOSAL: Use solo.qli.utils.log_time().
%
% PROPOSAL: Have local function db_get_ts() normalize data returned from
%           solo.db_get_ts() to TSeries, also for absent data.
%   NOTE: Would require argument for dimensions when empty.
%   PRO: Could (probably) simplify plot code a lot.
%
% PROPOSAL: Move quicklooks_24_6_2.m constants here. Submit values as arguments.
%
%
% quicklooks_24_6_2_h.m(), quicklooks_7day()
% ==========================================
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

% ##############################################################################



%=========================
% Assertions on arguments
%=========================
assert(islogical(generateNonweeklyQuicklooks) & islogical(generateNonweeklyQuicklooks))
assert(islogical(generateWeeklyQuicklooks   ) & islogical(generateWeeklyQuicklooks   ))
assert(iscolumn(DaysDtArray))
solo.qli.utils.assert_UTC_midnight_datetime(DaysDtArray)

% "Normalize"
DaysDtArray = unique(DaysDtArray);



%============
% ~Constants
%============
% IMPLEMENTATION NOTE: Disabling B (use empty; pretend there is no B data)
% speeds up solo.qli.quicklooks_24_6_2_h() greatly. Useful for some debugging.
% Should be enabled by default.
ENABLE_B                      = 1;    % 0 or 1.
% Whether to catch plotting exceptions, continue plotting other days/weeks, and
% then re-raise the last caught exception at the very end. This produces as many
% quicklooks as possible when one or some quicklooks fail.
% Should be enabled by default.
CATCH_PLOT_EXCEPTIONS_ENABLED = 1;    % 0 or 1.

% Path to IRF logo, relative to irfu-matlab root.
IRF_LOGO_RPATH = 'mission/solar_orbiter/+solo/irf_logo.png';

% NOTE: Usually found at /data/solo/data_yuri/.
VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
VHT_6H_DATA_FILENAME = 'V_RPW.mat';

% Define boundary of weeks. Beginning of stated weekday.
% ------------------------------------------------------
% NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
% Therefore using Wednesday as beginning of "week" for weekly quicklooks (until
% someone complains).
FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday



%==============================
% Miscellaneous initialization
%==============================
% Specify subdirectories for saving the respective types of quicklooks.
OutputPaths.path_2h  = fullfile(outputDir, '2h' );
OutputPaths.path_6h  = fullfile(outputDir, '6h' );
OutputPaths.path_24h = fullfile(outputDir, '24h');
OutputPaths.path_1w  = fullfile(outputDir, '1w' );

% Derive full path to IRF logo image if it seems that it should be used
% ---------------------------------------------------------------------
% Empty ==> Do not plot any logo.
% IMPLEMENTATION NOTE: The IRF logo should only be used for official quicklooks.
% The path is therefore only automatically set if it appears that the quicklooks
% being generated are "official".
irfLogoPath = [];
if isunix()
  [errorCode, stdoutStr] = system('hostname');
  assert(errorCode == 0, 'Error when calling "hostname". errorCode = %i', errorCode)
  hostName = strip(stdoutStr);

  if ismember(hostName, {'brain', 'spis'})
    parentDir          = fileparts(mfilename('fullpath'));
    irfumatlabRootPath = fullfile(parentDir, '../../../../');
    irfLogoPath        = fullfile(irfumatlabRootPath, IRF_LOGO_RPATH);
  end
end



tSec = tic();

% Array of plotting exceptions caught.
PlotExcArray = MException.empty(1, 0);

% Derive weeks from specified days (midnights which begin 7-day periods).
WeeksDtArray = solo.qli.utils.derive_weeks(DaysDtArray, FIRST_DAY_OF_WEEK);



%=========================================================
% Create subdirectories for different types of quicklooks
%=========================================================
% NOTE: Creates subdirectories for all four quicklook types, regardless of whether
%       they will actually be generated or not.
%   NOTE: This simplifies bash wrapper scripts that copy content of
%         sub-directories to analogous sub-directories, also when not all
%         quicklook types are generated.
for fnCa = fieldnames(OutputPaths)'
  dirPath = OutputPaths.(fnCa{1});
  [parentDir, dirBasename, dirSuffix] = fileparts(dirPath);
  % NOTE: Works (without warnings) also if subdirectories pre-exist ("msg"
  % contains warning which is never printed.)
  [success, msg] = mkdir(parentDir, [dirBasename, dirSuffix]);
  assert(success, 'Failed to create directory "%s". %s', dirPath, msg)
end



%=============================================
% Run the code for 2-, 6-, 24-hour quicklooks
%=============================================
if generateNonweeklyQuicklooks
  % Daily time-intervals

  % Load VHT data
  % -------------
  % This is the .mat file containing RPW speeds at 1h resolution.
  % The file should be in the current path. This file can be found in
  % brain:/solo/data/data_yuri/.
  vht1h = load(fullfile(vhtDataDir, VHT_1H_DATA_FILENAME));

  for iDay = 1:length(DaysDtArray)
    DayDt = DaysDtArray(iDay);

    try
      quicklooks_24_6_2_h_local(DayDt, vht1h, OutputPaths, irfLogoPath, ENABLE_B)
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(CATCH_PLOT_EXCEPTIONS_ENABLED, Exc)
    end
  end
end



%===================================
% Run the code for weekly overviews
%===================================
if generateWeeklyQuicklooks

  % Load VHT data
  % -------------
  % This is the .mat file containing RPW speeds at 6h resolution.
  % The file should be in the same folder as this script (quicklook_main).
  vht6h = load(fullfile(vhtDataDir, VHT_6H_DATA_FILENAME));

  for iWeek = 1:numel(WeeksDtArray)
    WeekDt = WeeksDtArray(iWeek);

    try
      quicklooks_7days_local(WeekDt, vht6h, OutputPaths, irfLogoPath)
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(CATCH_PLOT_EXCEPTIONS_ENABLED, Exc)
    end
  end
end



wallTimeSec     = toc(tSec);
wallTimeHours   = wallTimeSec/3600;
nQuicklooksDays = numel(DaysDtArray);

% NOTE: Execution speed may vary by orders of magnitude depending on settings
% (nonweekly vs weekly plots). May therefore want scientific notation.
fprintf('Wall time used:                       %g [h] = %g [s]\n', wallTimeHours, wallTimeSec);
fprintf('Wall time used per day of quicklooks: %g [h/day]\n',      wallTimeHours / nQuicklooksDays);



if CATCH_PLOT_EXCEPTIONS_ENABLED && ~isempty(PlotExcArray)
  fprintf(2, 'Caught %i plotting exceptions.\n', numel(PlotExcArray))
  fprintf(2, 'Rethrowing old (last) exception.\n')
  % NOTE: This does display (stderr) the stack trace for position
  % of the *ORIGINAL* error.
  rethrow(PlotExcArray(end))
end

end    % function



% Log time interval for which a plotting function is called.
%
% This is useful for more easily determining for which time interval the code
% crashes by reading the log.
function log_plot_function_time_interval(Tint)
utcStr1 = Tint(1).utc;
utcStr2 = Tint(2).utc;
% NOTE: Truncating subseconds (keeping accuracy down to seconds).
utcStr1 = utcStr1(1:19);
utcStr2 = utcStr2(1:19);

% Not specifying which plot function is called (weekly, nonweekly plots).
fprintf('Calling plot function for %s--%s.\n', utcStr1, utcStr2);
end



% Handle *PLOTTING* exception.
%
% Historically, the plotting code has caused many exceptions. One may want
% different behaviour depending on context.
%
% Production: Produce as many plots as possible.
%             => Catch exception and continue.
% Testing:    Crash on first exception so that it can be fixed.
%             => Rethrow exception as soon as possible.
function handle_plot_exception(catchExceptionEnabled, Exc)
if catchExceptionEnabled
  % Print stack trace without rethrowing exception.
  % One wants that in the log.
  % NOTE: fprintf(FID=2) => stderr
  fprintf(2, 'Caught plotting error without rethrowing it.\n')
  fprintf(2, 'Plot error/exception: "%s"\n', Exc.message)
  for i = 1:numel(Exc.stack)
    s = Exc.stack(i);
    fprintf(2, '    Error in %s (line %i)\n', s.name, s.line)
  end
else
  rethrow(Exc)
end
end



function quicklooks_24_6_2_h_local(Dt, vht1h, OutputPaths, irfLogoPath, enableB)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(1))
];
log_plot_function_time_interval(Tint)

Data = [];

Data.Vrpw   = vht1h.V_RPW_1h.tlim(Tint);
% E-field
Data.E      = db_get_ts(     'solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);
% RPW density
Data.Ne     = db_get_ts(     'solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);
% B-field
Data.B      = db_get_ts(     'solo_L2_mag-rtn-normal', 'B_RTN', Tint);
% Proton & alpha temperature
Data.Tpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity
Data.Vpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density
Data.Npas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux', 'eflux',Tint);
% TNR E-field
Data.Etnr   = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);
% Solar Orbiter position
% NOTE: Uses SPICE kernels indirectly. Kernels should be taken care of by
% "solo.get_position".
Data.solopos = get_SolO_position(Tint);

% Earth position (also uses SPICE)
DT = 60*60;
Data.earthpos = get_Earth_position(Tint, DT);

if ~enableB
  Data.B = [];
end

% Plot data and save figure
solo.qli.quicklooks_24_6_2_h(Data, OutputPaths, Tint, irfLogoPath)
end



function quicklooks_7days_local(Dt, vht6h, OutputPaths, irfLogoPath)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(7)), ...
];
log_plot_function_time_interval(Tint)

Data = [];

Data.Vrpw   = vht6h.V_RPW.tlim(Tint);
% E-field:
Data.E      = db_get_ts(     'solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);
% RPW density:
Data.Ne     = db_get_ts(     'solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);
% B-field:
Data.B      = db_get_ts(     'solo_L2_mag-rtn-normal-1-minute', 'B_RTN', Tint);
% Proton & alpha temperature:
Data.Tpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity:
Data.Vpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density:
Data.Npas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux', 'eflux',Tint);
% TNR E-field
Data.Etnr   = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);
% Solar Orbiter position
% NOTE: Uses SPICE kernels indirectly. Kernels should be taken care of by
% "solo.get_position()".
Data.solopos = get_SolO_position(Tint);

% Earth position (also uses SPICE)
DT = 60*60;
earthPosTSeries = get_Earth_position(Tint, DT);
Data.earthpos   = earthPosTSeries;

% Plot data and save figure
solo.qli.quicklooks_7days(Data, OutputPaths, Tint, irfLogoPath)
end



% Get Solar Orbiter position
%
% NOTE: Uses SPICE and "solo.get_position()" which can itself load SPICE
% kernels(!).
function soloPos = get_SolO_position(Tint)
assert((length(Tint) == 2) & isa(Tint, 'EpochTT'))

% IM = irfu-matlab (as opposed to SPICE).
imSoloPos = solo.get_position(Tint, 'frame', 'ECLIPJ2000');

% BUG?!!: If solo.get_position() is non-empty, and presumably contains a
%         value, THEN use SPICE value anyway?!! Note: This behaviour does
%         however mimick the behaviour of the original code by
%         Konrad Steinvall, IRF (before refactoring).
if ~isempty(imSoloPos)
  [radius, lon, lat] = cspice_reclat(imSoloPos.data');
  soloPos = irf.ts_vec_xyz(imSoloPos.time,[radius',lon',lat']);
else
  soloPos = imSoloPos;
end
end



% Use SPICE to get Earth's position.
function earthPosTSeries = get_Earth_position(Tint, dt)
assert(length(Tint) == 2)
assert(isnumeric(dt))

et = Tint.start.tts : dt : Tint.stop.tts;

spiceEarthPos = cspice_spkpos('Earth', et, 'ECLIPJ2000', 'LT+s', 'Sun');
if ~isempty(spiceEarthPos)
  [E_radius, E_lon, E_lat] = cspice_reclat(spiceEarthPos);
  earthPos = [E_radius', E_lon', E_lat'];

  Tlength  = Tint(end)-Tint(1);
  dTimes   = 0:dt:Tlength;
  Times    = Tint(1)+dTimes;
  earthPosTSeries = irf.ts_vec_xyz(Times, earthPos);
else
  % earthPos = [];
  earthPosTSeries = TSeries();   % Empty TSeries.
end
end



% Wrapper around solo.db_get_ts() which normalizes the output to always return
% one TSeries object.
%
% NOTE: solo.db_get_ts() returns a cell array of TSeries instead of a single
% TSeries when the underlying code thinks that the underlying CDFs do not have
% consistent metadata. See solo.db_get_ts().
%
function Ts = db_get_ts(varargin)

temp = solo.db_get_ts(varargin{:});

% Normalize (TSeries or cell array) --> TSeries.
if iscell(temp)
  temp = cell_array_TS_to_TS(temp);
end

Ts = temp;
end



% Take a cell array of TSeries and merges them into one TSeries.
function OutputTs = cell_array_TS_to_TS(InputTs)
assert(iscell(InputTs))

nCells   = numel(InputTs);
OutputTs = InputTs{1};

if nCells>1
  for iCell = 2:nCells    % NOTE: Begins at 2.
    OutputTs = OutputTs.combine(InputTs{iCell});
  end
end
end
