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
%   lookup before calling this function.
% * Uses SPICE implicitly, and therefore relies on some path convention used by
%   irfu-matlab for where to find SPICE kernels. Not sure which convention, but
%   presumably it does at least find /data/solo/SPICE/.
% * Uses solo.read_TNR() indirectly which in turns relies on a hardcoded
%   path to "/data/solo/remote/data/L2/thr/" and selected subdirectories.
% * Creates subdirectories to the output directory if not pre-existing.
% * Note: 7-day quicklooks always begin with a specific hardcoded weekday
%   (Wednesday as of 2024-03-07).
% * Overwrites pre-existing quicklook files without warning.
% * ~BUG: SolO DB (solo.db_get_ts() etc.) requires the caller to specify
%    dataset_ID plus "-cdag" if present, but does not raise exception if wrong.
%    ==> The code requires the datasets searched by SolO DB to have/lack "-cdag"
%    exactly as specified (hardcoded). If they are not, then data appears to not
%    be present and no exception is raised! This means that the code might not
%    recognize datasets for the path specified with solo.db_init().
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
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function generate_quicklooks_all_types(...
  irfLogoPath, vhtDataDir, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
%
% NOTE: Mission begins on 2020-02-12=Wednesday.
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
% PROPOSAL: Log time consumption for each call to plot functions.
%   PROPOSAL: Use solo.qli.utils.log_time().
%
% PROPOSAL: Parallelize loop over days and weeks.
%
% PROBLEM: How handle IRF logo file?!
%   NOTE: Is not in the git repo itself, only the directory structure.
%   PROPOSAL: Argument to generate_quicklooks* functions. cron code hardcodes
%             it.
%       TODO-DEC: Which default location to use (for cron.*)?
%           NOTE: Should work for both nas and locally.
%           PROPOSAL: In git repo. -- IMPLEMENTED
%               PRO: Connects file with s/w.
%               CON: File is not part of git repo.
%               PROPOSAL: Directly under git repo root.
%           PROPOSAL: Under home directory.
%               CON: Risks forgetting why it is needed to be there.
%               CON: No direct connection between QLI s/w and file.
%
% PROPOSAL: Have local function db_get_ts() normalize data returned from
%           solo.db_get_ts() to TSeries, also for absent data.
%   NOTE: Would require argument for dimensions when empty.
%   PRO: Could (probably) simplify plot code a lot.
%
% PROPOSAL: Always call solo.db_get_ts() both with and without "-cdag", to make
%           sure that the code does select non-existing datasets.
%           Use pre-existing db_get_ts() wrapper.
%
% PROPOSAL: Convert
%           generate_quicklooks_24h_6h_2h_local()
%           generate_quicklook_7days_local()
%           into separate function files which use SolO DB for retrieving values.
%     PRO: Usage of SolO DB becomes clearer for users.
%     PRO: It becomes easier for users to add new variables.
%     PRO: Can still use automated tests.
%     CON/PROBLEM: They both use wrapper db_get_ts() (function in this file).
%         CON-PROPOSAL: Move db_get_ts() to separate function and call it from
%                       above functions.
%         CON-PROPOSAL: Have
%                 solo.qli.generate_quicklooks_24h_6h_2h()
%                 solo.qli.generate_quicklook_7days()
%                 normalize their arguments with cell_array_TS_to_TS().
%             CON: Can not extend wrapper db_get_ts() to retrieve data for both
%                  -cdag and non-cdag.
%     PROPOSAL: Names *_local --> *_SolO_DB
%     CON/PROBLEM: They both use get_Earth_position(), get_SolO_position().
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
assert(islogical(generateNonweeklyQuicklooks) & isscalar(generateNonweeklyQuicklooks))
assert(islogical(generateWeeklyQuicklooks   ) & isscalar(generateWeeklyQuicklooks   ))
assert(iscolumn(DaysDtArray))
solo.qli.utils.assert_UTC_midnight_datetime(DaysDtArray)

% "Normalize"
DaysDtArray = unique(DaysDtArray);



%==============================
% Miscellaneous initialization
%==============================
OutputPaths = solo.qli.utils.create_output_directories(outputDir);

% Try to determine whether quicklooks are being generated as part of IRFU's
% official processing.
isOfficialProcessing = false;
if isunix()
  [errorCode, stdoutStr] = system('hostname');
  assert(errorCode == 0, 'Error when calling "hostname". errorCode = %i', errorCode)
  hostName = strip(stdoutStr);

  if ismember(hostName, solo.qli.const.OFFICIAL_PROCESSING_IRFU_HOST_NAMES_CA)
    isOfficialProcessing = true;
  end
end



tSec = tic();

% Array of plotting exceptions caught.
PlotExcArray = MException.empty(1, 0);

% Derive weeks from specified days (midnights which begin 7-day periods).
WeeksDtArray = solo.qli.utils.derive_weeks(DaysDtArray, solo.qli.const.FIRST_DAY_OF_WEEK);



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
  vht1h = load(fullfile(vhtDataDir, solo.qli.const.VHT_1H_DATA_FILENAME));

  for iDay = 1:length(DaysDtArray)
    DayDt = DaysDtArray(iDay);

    try
      trigger_automount(isOfficialProcessing)
      generate_quicklooks_24h_6h_2h_local(DayDt, vht1h, OutputPaths, irfLogoPath)
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

  % Load VHT data
  % -------------
  % This is the .mat file containing RPW speeds at 6h resolution.
  % The file should be in the same folder as this script (quicklook_main).
  vht6h = load(fullfile(vhtDataDir, solo.qli.const.VHT_6H_DATA_FILENAME));

  for iWeek = 1:numel(WeeksDtArray)
    WeekDt = WeeksDtArray(iWeek);

    try
      trigger_automount(isOfficialProcessing)
      generate_quicklook_7days_local(WeekDt, vht6h, OutputPaths, irfLogoPath)
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
fprintf('Wall time used:                       %g [h] = %g [s]\n', wallTimeHours, wallTimeSec);
fprintf('Wall time used per day of quicklooks: %g [h/day]\n',      wallTimeHours / nQuicklooksDays);



if solo.qli.const.CATCH_PLOT_EXCEPTIONS_ENABLED && ~isempty(PlotExcArray)
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
function handle_plot_exception(Exc)
if solo.qli.const.CATCH_PLOT_EXCEPTIONS_ENABLED
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



function generate_quicklooks_24h_6h_2h_local(Dt, vht1h, OutputPaths, irfLogoPath)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(1))
  ];
log_plot_function_time_interval(Tint)

Data = [];

Data.Vrpw   = vht1h.V_RPW_1h.tlim(Tint);
% E-field
Data.E      = db_get_ts(     'solo_L3_rpw-bia-efield-10-seconds-cdag', 'EDC_SRF', Tint);
% RPW density
Data.Ne     = db_get_ts(     'solo_L3_rpw-bia-density-10-seconds-cdag', 'DENSITY', Tint);
% B-field
Data.B      = db_get_ts(     'solo_L2_mag-rtn-normal', 'B_RTN', Tint);
% Proton & alpha temperature
Data.Tpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity
Data.Vpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density
Data.Npas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux', 'eflux', Tint);

% TNR E-field
% -----------
% BUG? Is not anything like an "E-field"!! Reading the wrong variable or
% mislabelling the right variable?
% NOTE: Variable is not used very much. Code only checks if empty or not.
%
%      FIELDNAM        (CDF_CHAR/8): "TNR_BAND"
%      CATDESC         (CDF_CHAR/31): "TNR band of the current record "
%      VAR_NOTES       (CDF_CHAR/71): "TNR band of the current record. Possible values are: 1=A, 2=B, 3=C, 4=D"
% /solo_L2_rpw-tnr-surv-cdag_20240101_V02.cdf
%
Data.Etnr   = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);
% Solar Orbiter position
% NOTE: Uses SPICE kernels indirectly. Kernels should be taken care of by
% solo.get_position().
Data.soloPos = get_SolO_position(Tint);

% Earth position (also uses SPICE)
DT = 60*60;
Data.earthPos = get_Earth_position(Tint, DT);

if ~solo.qli.const.ENABLE_B
  Data.B = [];
end

% Plot data and save figure
solo.qli.generate_quicklooks_24h_6h_2h(Data, OutputPaths, Tint, irfLogoPath)
end



function generate_quicklook_7days_local(Dt, vht6h, OutputPaths, irfLogoPath)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(7)), ...
  ];
log_plot_function_time_interval(Tint)

Data = [];

Data.Vrpw   = vht6h.V_RPW.tlim(Tint);
% E-field:
Data.E      = db_get_ts(     'solo_L3_rpw-bia-efield-10-seconds-cdag', 'EDC_SRF', Tint);
% RPW density:
Data.Ne     = db_get_ts(     'solo_L3_rpw-bia-density-10-seconds-cdag', 'DENSITY', Tint);
% B-field:
Data.B      = db_get_ts(     'solo_L2_mag-rtn-normal-1-minute', 'B_RTN', Tint);
% Proton & alpha temperature:
Data.Tpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity:
Data.Vpas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density:
Data.Npas   = db_get_ts(     'solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux', 'eflux', Tint);
% TNR E-field
Data.Etnr   = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);
% Solar Orbiter position
% NOTE: Uses SPICE kernels indirectly. Kernels should be taken care of by
% "solo.get_position()".
Data.soloPos = get_SolO_position(Tint);

% Earth position (also uses SPICE)
DT = 60*60;
earthPosTSeries = get_Earth_position(Tint, DT);
Data.earthPos   = earthPosTSeries;

% Plot data and save figure
solo.qli.generate_quicklook_7days(Data, OutputPaths, Tint, irfLogoPath)
end



% Use SPICE to get Solar Orbiter's position as
% [soloSunDistance, soloEclLongitude, soloEclLatitude].
%
% NOTE: Uses SPICE kernels and "solo.get_position()" indirectly through
% irfu-matlab which can itself load SPICE kernels(!).
function soloPosRadLonLatTSeries = get_SolO_position(Tint)
assert((length(Tint) == 2) & isa(Tint, 'EpochTT'))

% IM = irfu-matlab (as opposed to SPICE).
% See get_Earth_position() (in this file) for information on the coordinate
% system.
soloPosXyz = solo.get_position(Tint, 'frame', 'ECLIPJ2000');

if ~isempty(soloPosXyz)
  [soloSunDistance, soloEclLongitude, soloEclLatitude] = cspice_reclat(soloPosXyz.data');
  soloPosRadLonLatTSeries = irf.ts_vec_xyz(soloPosXyz.time, [soloSunDistance', soloEclLongitude', soloEclLatitude']);
else
  %soloPosRadLonLat = soloPosXyz;
  soloPosRadLonLatTSeries = TSeries();   % Empty TSeries.
end
end



% Use SPICE to get Earth's position as
% [earthSunDistance, earthEclLongitude, earthEclLatitude].
%
function earthPosRadLonLatTSeries = get_Earth_position(Tint, dt)
%===============================================================================
% Arguments for cspice_spkpos()
% -----------------------------
% 17  ECLIPJ2000  Ecliptic coordinates based upon the
%                 J2000 frame.
%
%                 The value for the obliquity of the
%                 ecliptic at J2000 is taken from page 114
%                 of [7] equation 3.222-1. This agrees with the
%                 expression given in [5].
%
% Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
% --
% 'LT+S'     Correct for one-way light time and
%            stellar aberration using a Newtonian
%            formulation. This option modifies the
%            position obtained with the 'LT' option
%            to account for the observer's velocity
%            relative to the solar system
%            barycenter. The result is the apparent
%            position of the target---the position
%            as seen by the observer.
%
% Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_spkpos.html
%===============================================================================
assert(length(Tint) == 2)
assert(isnumeric(dt))

et = Tint.start.tts : dt : Tint.stop.tts;

earthPosXyz = cspice_spkpos('Earth', et, 'ECLIPJ2000', 'LT+s', 'Sun');

if ~isempty(earthPosXyz)
  [earthSunDistance, earthEclLongitude, earthEclLatitude] = cspice_reclat(earthPosXyz);
  earthPos = [earthSunDistance', earthEclLongitude', earthEclLatitude'];

  Tlength = Tint(end)-Tint(1);
  dTimes  = 0:dt:Tlength;
  Times   = Tint(1)+dTimes;
  earthPosRadLonLatTSeries = irf.ts_vec_xyz(Times, earthPos);
else
  earthPosRadLonLatTSeries = TSeries();   % Empty TSeries.
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



% Try to trigger automount for official processing, to avoid file-reading
% problems for long runs.
%
% CONTEXT
% =======
% The automounting which should make /data/solo/ available on brain/spis (IRFU)
% does not always work as intended. It appears to not always work fast enough
% for reading files on rare occassions, leading to apparent can-not-read-file
% errors for existing files when making long runs, in particular on
% out-of-office hours(?). Other code which has faced this problem has
% seemingly resolved it by using Linux commands to list the directory contents
% and ignore the result before using /data/solo/ (such as
% "ls /data/solo/ >> /dev/null").
%
function trigger_automount(isOfficialProcessing)
if isOfficialProcessing
  junk = dir(solo.qli.const.OFFICIAL_PROCESSING_AUTOMOUNT_DIR);
end
% NOTE: Command only works on UNIX/Linux. Not Windows, MacOs etc.
% errorCode = system('ls -l /data/solo/ >> /dev/null');
end
