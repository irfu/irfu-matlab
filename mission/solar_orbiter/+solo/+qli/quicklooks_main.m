%
% Generate multiple IRFU-local types of quicklooks (files) for SolO data (not
% just RPW) for a specified time interval.
%
% Intended for batch processing, e.g. being called from bash script, e.g. cron
% job via a MATLAB wrapper script.
%
% NOTE: Requires solo.db_init() to have been properly used to initialize dataset
%       lookup.
% NOTE: Uses SPICE implicitly, and therefore relies on some path convention. Not
%       sure which, but presumably it does at least find /data/solo/SPICE/.
% NOTE: Uses solo.read_TNR() indirectly which in turns relies on a hardcoded
%       path to "/data/solo/remote/data/L2/thr/" and selected subdirectories.
% NOTE: Creates subdirectories to the output directory if not pre-existing.
% NOTE: 7-day plots will only cover that largest sub-time interval that is
%       composed of full 7-day periods. Note: 7-day periods begin with a
%       hardcoded weekday (Wednesday as of 2023-05-09).
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
%       Whether to run the resp. groups of plots.
%       NOTE: Permits chars "0" and "1" for when calling from bash.
%       Useful for testing and not re-running unnecessary time-consuming plots.
% utcBegin, utcEnd : Strings.
%       Defines time interval for which quicklooks should be generated.
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson, IRF, Uppsala, Sweden.
%
 function quicklooks_main(...
        logoPath, vhtDataDir, outputDir, ...
        runNonweeklyPlots, runWeeklyPlots, utcBeginStr, utcEndStr)
%
% NOTE: Mission begins on 2020-02-12=Wednesday.
% ==> There is no SPICE data on Monday-Tuesday before this date.
% ==> Code fails for week Monday-to-Sunday.
%     PROPOSAL: Additionally round start time up to start date.
%
% TODO-NI: Overwrites pre-existing plots?
%
% PROPOSAL: Better name: quicklooks_main?
%       ~main, ~plot, ~generate
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
% PROPOSAL: Catch exceptions raised by plotting subfunctions?
%   PROPOSAL: Make configurable.
%   PROPOSAL: Config 1: Catch exceptions but return error code after completing
%             entire time interval.
%       NOTE: Can catch multiple crashes/bugs in one run.
%       PRO: Good for debugging since can find multiple dates with crashes in
%            one (long, e.g. overnight) test run.
%       PRO: Useful for batch processing, cron jobs.
%           PRO: Not worse than the old crash immediately behaviour.
%   PROPOSAL: Config 2: Immediately raise exception when catching exception.
%       PRO: Good for debugging (?).
%   PROPOSAL: Config 3: Catch all exceptions.
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
% PROPOSAL: Log time period being plotted.
%   PRO: Useful for identifying time period/day that causes crash.
% PROPOSAL: Log time consumption for each call to plot functions.
%   PROPOSAL: Use solo.qli.utils.log_time().
%
%
% quicklooks_24_6_2_h.m(), quicklooks_7day()
% ==========================================
% PROBLEM: Lots of cases of, and checks for data that may or may not be missing.
% PROBLEM: Missing data can be represented as [] (0x0), rather than e.g. Nx0.
%   Stems from solo.db_get_ts() (?) not finding any data for selected time
%   interval. Can not know the dimensions of missing data if can not find any
%   CDF at all.
% PROPOSAL: Function for normalizing data returned from solo.db_get_ts() to
%           zero-length TSeries.
%   CON: Can not do since does not know non-time dimensions.
% PROPOSAL: Iterations over 6h an 2h seem to be identical. Wrap in function(s).
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


%============
% ~Constants
%============
% IMPLEMENTATION NOTE: Disabling B (use empty; pretend there is no B data)
% speeds up solo.qli.quicklooks_24_6_2_h() greatly. Useful for some debugging.
ENABLE_B = 1;

% Specify subdirectories for saving the respective types of plots.
PATHS.path_2h  = fullfile(outputDir, '2h' );
PATHS.path_6h  = fullfile(outputDir, '6h' );
PATHS.path_24h = fullfile(outputDir, '24h');
PATHS.path_1w  = fullfile(outputDir, '1w' );

% NOTE: Usually found on solo/data_yuri.
VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
VHT_6H_DATA_FILENAME = 'V_RPW.mat';

% Define boundary of weeks. Beginning of stated weekday.
% ------------------------------------------------------
% NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
% Therefore using Wednesday as beginning of "week" for weekly plots (until
% someone complains).
FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday



tSec = tic();



%=================================
% Specify time interval for plots
%=================================
% Non-weekly plots.
TimeIntervalNonWeeks = irf.tint(utcBeginStr, utcEndStr);
% Weekly plots: Indirectly sets weekly boundaries.
TimeIntervalWeeks = derive_TimeIntervalWeeks(...
    TimeIntervalNonWeeks(1), TimeIntervalNonWeeks(2), FIRST_DAY_OF_WEEK);



%=======================
% Create subdirectories
%=======================
% NOTE: Creates subdirectories for all four plot types, regardless of whether
%       they will actually be generated or not.
%   NOTE: This simplifies bash wrapper scripts that copy content of
%         sub-directories to analogous sub-directories, also when not all
%         plot types are generated.
for fnCa = fieldnames(PATHS)'
    dirPath = PATHS.(fnCa{1});
    [parentDir, dirBasename, dirSuffix] = fileparts(dirPath);
    % NOTE: Works (without warnings) also if subdirectories pre-exist ("msg"
    % contains warning which is never printed.)
    [success, msg] = mkdir(parentDir, [dirBasename, dirSuffix]);
    assert(success, 'Failed to create directory "%s". %s', dirPath, msg)
end



%=============================================
% Run the code for 2-, 6-, 24-hour quicklooks
%=============================================
if runNonweeklyPlots

    Time1DayStepsArray = make_time_array(TimeIntervalNonWeeks, 1); % Daily time-intervals

    % Load data
    % This is the .mat file containing RPW speeds at 1h resolution.
    % The file should be in the current path. This file can be found in
    % brain:/solo/data/data_yuri/.
    vht1h = load(fullfile(vhtDataDir, VHT_1H_DATA_FILENAME));

    for iTint=1:length(Time1DayStepsArray)-1
        % Select time interval.
        Tint=irf.tint(Time1DayStepsArray(iTint), Time1DayStepsArray(iTint+1));

        quicklooks_24_6_2_h_local(Tint, vht1h, PATHS, logoPath, ENABLE_B)
    end
end



%===================================
% Run the code for weekly overviews
%===================================
if runWeeklyPlots

    Time7DayStepsArray = make_time_array(TimeIntervalWeeks, 7);% weekly time-intervals

    % Load data
    % This is the .mat file containing RPW speeds at 6h resolution.
    % The file should be in the same folder as this script (quicklook_main).
    vht6h = load(fullfile(vhtDataDir, VHT_6H_DATA_FILENAME));

    for iTint=1:length(Time7DayStepsArray)-1
        % Select time interval.
        Tint = irf.tint(Time7DayStepsArray(iTint), Time7DayStepsArray(iTint+1));

        quicklooks_7days_local(Tint, vht6h, PATHS, logoPath)
    end
end



wallTimeSec   = toc(tSec);
wallTimeHours = wallTimeSec/3600;
plotsTimeDays = (TimeIntervalNonWeeks.tts(2) - TimeIntervalNonWeeks.tts(1)) / 86400;

% NOTE: Execution speed may vary by orders of magnitude depending on settings
% (nonweekly vs weekly plots). May therefore want scientific notation.
fprintf('Wall time used:                  %g [h] = %g [s]\n', wallTimeHours, wallTimeSec);
fprintf('Wall time used per day of plots: %g [h/day]\n',      wallTimeHours / plotsTimeDays);

end    % function



function quicklooks_24_6_2_h_local(Tint, vht1h, Paths, logoPath, enableB)
    Data = [];

    Data.Vrpw = vht1h.V_RPW_1h.tlim(Tint);

    % E-field:
    Data.E = db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    Data.Ne = db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    Data.B = db_get_ts('solo_L2_mag-rtn-normal','B_RTN', Tint);

    % Proton & alpha temperature:
    Data.Tpas = db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    Data.Vpas = db_get_ts('solo_L2_swa-pas-grnd-mom','V_RTN', Tint);

    % Proton & alpha density:
    Data.Npas = db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

    % Ion spectrum
    Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux','eflux',Tint);

    % TNR E-field
    Data.Etnr = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);

    % Solar Orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".
    Data.solopos = get_SolO_pos(Tint);

    % Earth position (also uses SPICE)
    dt = 60*60;
    Data.earthpos = get_Earth_pos(Tint, dt);

    if ~enableB
        Data.B = [];
    end

    % Plot data and save figure
    solo.qli.quicklooks_24_6_2_h(Data, Paths, Tint, logoPath)
end



function quicklooks_7days_local(Tint, vht6h, Paths, logoPath)
    Data = [];

    Data.Vrpw = vht6h.V_RPW.tlim(Tint);

    % E-field:
    Data.E = db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    Data.Ne = db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    Data.B = db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);

    % Proton & alpha temperature:
    Data.Tpas = db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    Data.Vpas = db_get_ts('solo_L2_swa-pas-grnd-mom','V_RTN', Tint);

    % Proton & alpha density:
    Data.Npas = db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

    % Ion spectrum
    Data.ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux','eflux',Tint);

    % TNR E-field
    Data.Etnr = solo.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);

    % Solar Orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".
    Data.solopos = get_SolO_pos(Tint);

    % Earth position (also uses SPICE)
    dt       = 60*60;
%     Tlength  = Tint(end)-Tint(1);
%     dTimes   = 0:dt:Tlength;
%     Times    = Tint(1)+dTimes;
    earthPosTSeries = get_Earth_pos(Tint, dt);
%     if ~isempty(earthPos)
%         Data.earthpos = irf.ts_vec_xyz(Times, earthPos);
%     else
%         Data.earthpos = TSeries();
%     end
    Data.earthpos = earthPosTSeries;

    % Plot data and save figure
    solo.qli.quicklooks_7days(Data, Paths, Tint, logoPath)
end



% Get Solar Orbiter position
%
% NOTE: Uses SPICE and "solo.get_position()".
function soloPos = get_SolO_pos(Tint)
    assert(length(Tint) == 2)
    
    % IM = irfu-matlab (as opposed to SPICE).
    imSoloPos = solo.get_position(Tint,'frame','ECLIPJ2000');

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
function earthPosTSeries = get_Earth_pos(Tint, dt)
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
        earthPosTSeries = TSeries();
    end
end



% Function for deriving the exact week boundaries to use.
function TimeIntervalWeeks = derive_TimeIntervalWeeks(TimeBegin, TimeEnd, firstDayOfWeek)
    assert(isscalar(TimeBegin))
    assert(isscalar(TimeEnd))

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
% TintInterval
%       Time interval
% nDays
%       Number of days between each timestamp.
function TimeArray = make_time_array(TintInterval, nDays)
    assert(length(TintInterval) == 2)

    t0          = TintInterval(1);
    tlength     = TintInterval(2) - TintInterval(1);
    % NOTE: Does not take leap seconds into account.
    stepSizeSec = nDays*24*60*60;    % seconds.

    dt          = 0:stepSizeSec:tlength;
    TimeArray   = t0+dt;
end



% Wrapper around solo.db_get_ts() that normalizes the output to a TSeries.
%
% NOTE: solo.db_get_ts() has been observed to return cell array of TSeries
% (instead of TSeries) for Npas, Tpas and Vpas for
% solo.qli.quicklooks_main('2020-10-21T00:00:00', '2020-10-28T00:00:00', '/data/solo/data_yuri', ...)
% There might be other cases but those are as of yet unknown.
% /Erik P G Johansson 2021-03-22
%
function Ts = db_get_ts(varargin)
    temp = solo.db_get_ts(varargin{:});

    % Normalize (TSeries or cell array) --> TSeries.
    if iscell(temp)
        temp = cell_array_TS_to_TS(temp);
    end

    Ts = temp;
end



% Takes a cell-array of TSeries and merges them to one TSeries.
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
% t1, t2 : GenericTimeArray, scalar.
function t2 = round_to_week(t1, roundDir, firstDayOfWeek)
    assert(isscalar(t1))
    assert(ismember(roundDir, [-1, 1]))



    dv1  = irf.cdf.TT2000_to_datevec(t1.ttns);
    dt1a = datetime(dv1, 'TimeZone', 'UTCLeapSeconds');

    % Round to midnight.
    dt1b = dateshift(dt1a, 'start', 'day');
    if (roundDir == 1) && (dt1a ~= dt1b)
        % IMPLEMENTATION NOTE: dateshift(..., 'end', 'day') "rounds" to one day after if
        % timestamp is already midnight. Therefore do not want use that.
        dt1b = dt1b + days(1);
    end

    % Round to week boundary, as defined by beginningOfWeek.
    % NOTE: dateshift( 'dayofweek' ) rounds to next match, including potentially
    % the same day.
    dt1c = dateshift(dt1b, 'dayofweek', firstDayOfWeek);
    if (roundDir == -1) && (dt1b ~= dt1c)
        dt1c = dt1c - days(7);
    end

    tt2000 = irf.cdf.datevec_to_TT2000(datevec(dt1c));
    t2     = irf.time_array(tt2000);
end
