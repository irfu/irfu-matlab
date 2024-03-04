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
% * 7-day plots will only cover that largest sub-time interval that is
%   composed of full 7-day periods. Note: 7-day periods begin with a
%   hardcoded weekday (Wednesday as of 2023-07-24).
% * Overwrites pre-existing plot files without warning.
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
%       Path to directory containing VHT (velocity) .mat files
%       V_RPW_1h.mat and V_RPW.mat. Typically brain:/data/solo/data_yuri/
% outputDir
%       Plots will be placed in subdirectories under this directory.
%       NOTE: Will create subdirectories if not pre-existing.
% generateNonweeklyPlots, generateWeeklyPlots
%       Whether to generate non-weekly (2h, 6h, 24h) quicklooks and or weekly
%       quicklooks.
%       Useful for testing and not re-running unnecessary time-consuming plots.
% utcBegin, utcEnd
%       Strings. Defines time interval for which quicklooks should be generated.
%       NOTE: Weekly plots will only be produced for those weeks which are
%       contained entirely inside the specified time interval (?). Weekly plots
%       always (as of 2023-07-24) begin on a Wednesday(!) at 00:00:00.
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function quicklooks_main(...
  logoPath, vhtDataDir, outputDir, ...
  generateNonweeklyPlots, generateWeeklyPlots, utcBeginStr, utcEndStr)
%
% NOTE: Mission begins on 2020-02-12=Wednesday.
% ==> There is no SPICE data on Monday-Tuesday before this date.
% ==> Code fails for week Monday-to-Sunday.
%     PROPOSAL: Additionally round start time up to start date.
%
% PROPOSAL: Better name: quicklooks_main?
%       ~main, ~plot, ~generate, ~qli
%       generate_quicklooks()
%
% PROPOSAL: Directly generate arrays of timestamps for iterating over, instead
%           of via TimeIntervalNonWeeks and TimeIntervalWeeks.
% PROPOSAL: Print time intervals for which
%           solo.qli.quicklooks_24_6_2_h() and
%           solo.qli.quicklooks_7days() are called.
%   PRO: Text will be stored in logs (not created by this code, but by bash
%        wrapper scripts).
%   PRO: Useful for more easily determining for which time intervals the code
%        (those two functions) crashes.
%
% PROPOSAL: Some way of handling disk access error?
%   PROPOSAL: try-catch plot code once (weekly or non-weekly plot function).
%             Then try without catch a second time, maybe after delay.
%             If the first call fails due to disk access error, it might still
%             trigger automount which makes the second attempt succeed.
%
% PROPOSAL: Refactor some functions into separate function files with test code.
%   Ex: derive_TimeIntervalWeeks()
%       make_time_array()
%       round_to_week()
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
% PROPOSAL: Always round time up to full weeks for 7-day plots.
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
assert(islogical(generateNonweeklyPlots) & islogical(generateNonweeklyPlots))
assert(islogical(generateWeeklyPlots   ) & islogical(generateWeeklyPlots   ))



%============
% ~Constants
%============
% IMPLEMENTATION NOTE: Disabling B (use empty; pretend there is no B data)
% speeds up solo.qli.quicklooks_24_6_2_h() greatly. Useful for some debugging.
% Should be enabled by default.
ENABLE_B                      = 1;    % 0 or 1.
% Whether to catch plotting exceptions, continue plotting other days/weeks, and
% then re-raise the last caught exception at the very end. This produces as many
% plots as possible when one or some plots fail.
% Should be enabled by default.
CATCH_PLOT_EXCEPTIONS_ENABLED = 1;    % 0 or 1.



% NOTE: Usually found at /data/solo/data_yuri/.
VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
VHT_6H_DATA_FILENAME = 'V_RPW.mat';

% Define boundary of weeks. Beginning of stated weekday.
% ------------------------------------------------------
% NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
% Therefore using Wednesday as beginning of "week" for weekly plots (until
% someone complains).
FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday



%==============================
% Miscellaneous initialization
%==============================
% Specify subdirectories for saving the respective types of plots.
OutputPaths.path_2h  = fullfile(outputDir, '2h' );
OutputPaths.path_6h  = fullfile(outputDir, '6h' );
OutputPaths.path_24h = fullfile(outputDir, '24h');
OutputPaths.path_1w  = fullfile(outputDir, '1w' );



tSec = tic();

% Array of plotting exceptions caught.
PlotExcArray = MException.empty(1, 0);



%=================================
% Specify time interval for plots
%=================================
% Non-weekly plots.
TimeIntervalNonWeeks = irf.tint(utcBeginStr, utcEndStr);
% Weekly plots: Indirectly sets weekly boundaries.
TimeIntervalWeeks = derive_TimeIntervalWeeks(...
  TimeIntervalNonWeeks(1), TimeIntervalNonWeeks(2), FIRST_DAY_OF_WEEK);



%=========================================================
% Create subdirectories for different types of quicklooks
%=========================================================
% NOTE: Creates subdirectories for all four plot types, regardless of whether
%       they will actually be generated or not.
%   NOTE: This simplifies bash wrapper scripts that copy content of
%         sub-directories to analogous sub-directories, also when not all
%         plot types are generated.
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
if generateNonweeklyPlots
  % Daily time-intervals
  Time1DayStepsArray = make_time_array(TimeIntervalNonWeeks, 1);

  % Load data
  % This is the .mat file containing RPW speeds at 1h resolution.
  % The file should be in the current path. This file can be found in
  % brain:/solo/data/data_yuri/.
  vht1h = load(fullfile(vhtDataDir, VHT_1H_DATA_FILENAME));

  for iTint=1:length(Time1DayStepsArray)-1
    % Select time interval.
    Tint = irf.tint(Time1DayStepsArray(iTint), Time1DayStepsArray(iTint+1));
    try
      quicklooks_24_6_2_h_local(Tint, vht1h, OutputPaths, logoPath, ENABLE_B)
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(CATCH_PLOT_EXCEPTIONS_ENABLED, Exc)
    end
  end
end



%===================================
% Run the code for weekly overviews
%===================================
if generateWeeklyPlots

  % Weekly time-intervals
  Time7DayStepsArray = make_time_array(TimeIntervalWeeks, 7);

  % Load data
  % This is the .mat file containing RPW speeds at 6h resolution.
  % The file should be in the same folder as this script (quicklook_main).
  vht6h = load(fullfile(vhtDataDir, VHT_6H_DATA_FILENAME));

  for iTint=1:length(Time7DayStepsArray)-1
    % Select time interval.
    Tint = irf.tint(Time7DayStepsArray(iTint), Time7DayStepsArray(iTint+1));
    try
      quicklooks_7days_local(Tint, vht6h, OutputPaths, logoPath)
    catch Exc
      PlotExcArray(end+1) = Exc;
      handle_plot_exception(CATCH_PLOT_EXCEPTIONS_ENABLED, Exc)
    end
  end
end



wallTimeSec   = toc(tSec);
wallTimeHours = wallTimeSec/3600;
plotsTimeDays = (TimeIntervalNonWeeks.tts(2) - TimeIntervalNonWeeks.tts(1)) / 86400;

% NOTE: Execution speed may vary by orders of magnitude depending on settings
% (nonweekly vs weekly plots). May therefore want scientific notation.
fprintf('Wall time used:                  %g [h] = %g [s]\n', wallTimeHours, wallTimeSec);
fprintf('Wall time used per day of plots: %g [h/day]\n',      wallTimeHours / plotsTimeDays);



if CATCH_PLOT_EXCEPTIONS_ENABLED && ~isempty(PlotExcArray)
  fprintf(2, 'Caught %i plotting exceptions.\n', numel(PlotExcArray))
  fprintf(2, 'Rethrowing old (last) exception.\n')
  % NOTE: This does display (stderr) the stack trace for position
  % of the *ORIGINAL* error.
  rethrow(PlotExcArray(end))
end

end    % function



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



function quicklooks_24_6_2_h_local(Tint, vht1h, OutputPaths, logoPath, enableB)
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
solo.qli.quicklooks_24_6_2_h(Data, OutputPaths, Tint, logoPath)
end



function quicklooks_7days_local(Tint, vht6h, OutputPaths, logoPath)
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
% "solo.get_position".
Data.solopos = get_SolO_position(Tint);

% Earth position (also uses SPICE)
DT = 60*60;
earthPosTSeries = get_Earth_position(Tint, DT);
Data.earthpos = earthPosTSeries;

% Plot data and save figure
solo.qli.quicklooks_7days(Data, OutputPaths, Tint, logoPath)
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



% Function for deriving the exact week boundaries to use.
function TimeIntervalWeeks = derive_TimeIntervalWeeks(TimeBegin, TimeEnd, firstDayOfWeek)
assert(isscalar(TimeBegin) & isa(TimeBegin, 'EpochTT'))
assert(isscalar(TimeEnd)   & isa(TimeBegin, 'EpochTT'))

tWeeksBegin = round_to_week(TimeBegin,  1, firstDayOfWeek);
tWeeksEnd   = round_to_week(TimeEnd,   -1, firstDayOfWeek);

if tWeeksBegin <= tWeeksEnd
  TimeIntervalWeeks = irf.tint(tWeeksBegin, tWeeksEnd);
else
  % Empty week. ~Hackish.
  TimeIntervalWeeks = irf.tint(tWeeksBegin, tWeeksBegin);
end
end



% Generate array of timestamps with specific and constant frequency.
%
% ARGUMENTS
% =========
% TintInterval
%       Time interval
% nDays
%       Number of days between each timestamp.
%
function TimeArray = make_time_array(TintInterval, nDays)
assert((length(TintInterval) == 2) & isa(TintInterval, 'EpochTT'))

t0          = TintInterval(1);
tlength     = TintInterval(2) - TintInterval(1);
% NOTE: Does not take leap seconds into account.
stepSizeSec = nDays*24*60*60;    % seconds.

dt          = 0:stepSizeSec:tlength;
TimeArray   = t0+dt;
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



% Round timestamp down/up to beginning of week.
%
% ARGUMENTS
% =========
%     t1
%         Scalar EpochTT.
%     roundDir
%         Scalar number. -1 or 1.
%
function t2 = round_to_week(t1, roundDir, firstDayOfWeek)
assert(isscalar(t1) & isa(t1, 'EpochTT'))
assert(ismember(roundDir, [-1, 1]))

dv1  = irf.cdf.TT2000_to_datevec(t1.ttns);
Dt1a = datetime(dv1, 'TimeZone', 'UTCLeapSeconds');

% Round to midnight.
Dt1b = dateshift(Dt1a, 'start', 'day');
if (roundDir == 1) && (Dt1a ~= Dt1b)
  % IMPLEMENTATION NOTE: dateshift(..., 'end', 'day') "rounds" to one day
  % after if timestamp is already midnight. Therefore do not want use
  % that.
  % 2024-03-04: This if statement might be unnecessary.
  Dt1b = Dt1b + days(1);
end

% Round to week boundary, as defined by beginningOfWeek.
% NOTE: dateshift( 'dayofweek' ) rounds to next match, including potentially
% the same day.
Dt1c = dateshift(Dt1b, 'dayofweek', firstDayOfWeek);
if (roundDir == -1) && (Dt1b ~= Dt1c)
  Dt1c = Dt1c - days(7);
end

tt2000 = irf.cdf.datevec_to_TT2000(datevec(Dt1c));
t2     = irf.time_array(tt2000);
end
