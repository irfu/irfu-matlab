%
% Generate multiple IRFU-local types of quicklooks (files) for SolO data (not
% just RPW) for a specified time interval.
%
% Intended for batch processing, e.g. being called from bash script.
%
% NOTE: Uses SPICE implicitly, and therefore relies on some path convention. Not
% sure which, but presumably it does at least find /data/solo/SPICE/.
% NOTE: Still uses some hard-coded paths for locating datasets (solo.db_init()).
% NOTE: Creates subdirectories to the output directory if not pre-existing.
% NOTE: The time coverage of weekly plots uses 7-day periods that always start
% with utcBegin.
% NOTE: There will only be plots that cover full 7-day periods. If the global
% time interval is not a multiple of 7 days, then end of that global time
% interval will not be covered by a 7-day plot.
%
%
% ARGUMENTS
% =========
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
% Uppsala, Sweden. Modified by Erik P G Johansson.
%
function quicklook_main(vhtDataDir, outputDir, runNonweeklyPlots, runWeeklyPlots, utcBegin, utcEnd)
%
% PROPOSAL: Log wall time used.
%   NOTE: ~Can not time per plot.
%   PROPOSAL: Log wall time per day of data.
%
% TODO-NI: Overwrites pre-existing plots?
%
% PROPOSAL: Reorg to separate internal functions for non-weekly and weekly plots
% respectively.
%   NOTE: Would need to have arguments for debugging constants like DISABLE_B etc.
%
% PROPOSAL: Always start weeks with same weekday.
%   NOTE: Currently only plots weeks which are entirely covered by input time
%   interval. Weeks always start with start date.
%   PROBLEM: If rounding to nearest previous/successive Monday, then one will
%   always plot too much, or too little.
%       PROPOSAL: Argument/setting for policy.
%           ~all     : Always plot time interval. Round to larger time interval.
%                      Incomplete plots if global time interval does not fall on Monday 00:00.
%           ~smaller : Round to smaller time interval.
%       PROPOSAL: Always round to smaller time interval. The caller has to be
%       aware of time intervals if want to cover specifik weeks.
%   PROBLEM: Mission begins on 2020-02-12=Wednesday.
%   ==> There is no SPICE data on Monday-Tuesday before this date.
%   ==> Code fails for week Monday-to-Sunday.
%       PROPOSAL: Additionally round start time up to start date.


runNonweeklyPlots = interpret_argument_flag(runNonweeklyPlots);
runWeeklyPlots    = interpret_argument_flag(runWeeklyPlots);



% IMPLEMENTATION NOTE: Needed to make "DB" work. Necessary when calling from
% bash.
irf



%================================
% Specify Solar Orbiter database
%================================
% For data on the IRFU server the following two are sufficient. You also need
% V_RPW.mat and V_RPW_1h.mat, found on solo/data_yuri.
solo.db_init('local_file_db', '/data/solo/');
solo.db_init('local_file_db', '/data/solo/data_irfu');

% Setup cache
solo.db_init('db_cache_size_max', 4096)
solo.db_cache('on', 'save')



%===========
% Constants
%===========
% "Disabling B" speeds up solo.quicklooks_24_6_2_h(). Useful for testing.
ENABLE_B = 1;

% Specify subdirectories for saving the respective types of plots.
PATHS.path_2h  = fullfile(outputDir, '2h' );
PATHS.path_6h  = fullfile(outputDir, '6h' );
PATHS.path_24h = fullfile(outputDir, '24h');
PATHS.path_1w  = fullfile(outputDir, '1w' );

VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
VHT_6H_DATA_FILENAME = 'V_RPW.mat';

% Beginning of SPICE kernels.
% MISSION_BEGIN = irf.time_array('2020-02-12T00:00:00');



tSec = tic();

% Specify time interval for plots
TimeInterval      = irf.tint(utcBegin, utcEnd);
% TimeIntervalWeeks = irf.tint(...
%     round_to_week(TimeInterval(1),  1), ...
%     round_to_week(TimeInterval(2), -1));
TimeIntervalWeeks      = TimeInterval;


%=======================
% Create subdirectories
%=======================
for fnCa = fieldnames(PATHS)'
    dirPath = PATHS.(fnCa{1});
    [parentDir, dirBasename, dirSuffix] = fileparts(dirPath);
    % NOTE: Works (without warnings) also if subdirectories pre-exist ("msg"
    % contains warning which is enver printed.)
    [success, msg] = mkdir(parentDir, [dirBasename, dirSuffix]);
    assert(success, 'Failed to create directory "%s". %s', dirPath, msg)
end



%=============================================
% Run the code for 2-, 6-, 24-hour quicklooks
%=============================================
if runNonweeklyPlots
    
    times_1d = make_tints(TimeInterval, 1); % Daily time-intervals

    % Load data
    % This is the .mat file containing RPW speeds at 1h resolution.
    % The file should be in the current path. This file can be found in
    % /solo/data/data_yuri/.
    S = load(fullfile(vhtDataDir, VHT_1H_DATA_FILENAME));

    for iTint=1:length(times_1d)-1
        % Time interval
        Tint=irf.tint(times_1d(iTint), times_1d(iTint+1));
        
        Data.Vrpw = S.V_RPW_1h.tlim(Tint);
        %%
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
        
        % Solar Orbiter position
        % Note: solopos uses SPICE, but should be taken care of by
        % "solo.get_position".
        posSolO = solo.get_position(Tint,'frame','ECLIPJ2000');        
        if ~isempty(posSolO)
            [radius, lon, lat] = cspice_reclat(posSolO.data');
            Data.solopos = irf.ts_vec_xyz(posSolO.time,[radius',lon',lat']);
        else
            Data.solopos=posSolO;
        end
        
        % Earth position (also uses SPICE)
        dt=60*60;
        et = Tint.start.tts:dt:Tint.stop.tts;
        posEarth= cspice_spkpos('Earth', et, 'ECLIPJ2000', 'LT+s', 'Sun');
        if ~isempty(posEarth)
            [E_radius, E_lon, E_lat] = cspice_reclat(posEarth);
            Data.earthpos = [E_radius',E_lon',E_lat'];    
        else
            Data.earthpos=[];
        end
        
        
        if ~ENABLE_B
            Data.B = [];
        end

        % Plot data and save figure
        solo.quicklooks_24_6_2_h(Data,PATHS,Tint)
        
    end    % for
end



%===================================
% Run the code for weekly overviews
%===================================
if runWeeklyPlots
        
    times_7d = make_tints(TimeIntervalWeeks, 7);% weekly time-intervals
    
    % Load data
    % This is the .mat file containing RPW speeds at 6h resolution.
    % The file should be in the same folder as this script (quicklook_main).
    S = load(fullfile(vhtDataDir, VHT_6H_DATA_FILENAME));
    
    for iTint=1:length(times_7d)-1
  
        % Time interval
        Tint = irf.tint(times_7d(iTint), times_7d(iTint+1));
        
        Data2.Vrpw = S.V_RPW.tlim(Tint);
        
        % E-field:
        Data2.E = db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);
        
        % RPW density:
        Data2.Ne = db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);
        
        % B-field:
        Data2.B = db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);
        
        % Proton & alpha temperature:
        Data2.Tpas = db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);
        
        % Proton & alpha velocity:
        Data2.Vpas = db_get_ts('solo_L2_swa-pas-grnd-mom','V_RTN', Tint);
        
        % Proton & alpha density:
        Data2.Npas = db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);
        
        % Solar Orbiter position
        % Note: solopos uses SPICE, but should be taken care of by
        % "solo.get_position".
        posSolO = solo.get_position(Tint,'frame','ECLIPJ2000');        
        if ~isempty(posSolO)
            [radius, lon, lat] = cspice_reclat(posSolO.data');
            Data2.solopos = irf.ts_vec_xyz(posSolO.time,[radius',lon',lat']);
        else
            Data2.solopos=posSolO;
        end
        
        % Earth position (also uses SPICE)
        dt=60*60;
        et = Tint.start.tts:dt:Tint.stop.tts;
        Tlength=Tint(end)-Tint(1);
        dTimes = 0:dt:Tlength;
        Times = Tint(1)+dTimes;
        posEarth= cspice_spkpos('Earth', et, 'ECLIPJ2000', 'LT+s', 'Sun');
        if ~isempty(posEarth)
            [E_radius, E_lon, E_lat] = cspice_reclat(posEarth);
            Data2.earthpos = irf.ts_vec_xyz(Times,[E_radius',E_lon',E_lat']);    
        else
            Data2.earthpos=TSeries();
        end
        
        % Plot data and save figure
        solo.quicklooks_7days(Data2,PATHS,Tint)
        
    end    % for
end



wallTimeSec   = toc(tSec);
wallTimeHours = wallTimeSec/3600;
plotsTimeDays = (TimeInterval.tts(2) - TimeInterval.tts(1)) / 86400;

% NOTE: Execution speed may vary by orders of magnitude depending on settings
% (nonweekly vs weekly plots). May therefore want scientific notation.
fprintf('Wall time used:                  %g [h] = %g [s]\n',     wallTimeHours, wallTimeSec);
fprintf('Wall time used per day of plots: %g [h/day]\n', wallTimeHours / plotsTimeDays);

end    % function



% Auxilliary function
%
function out_times = make_tints(Tint, nDays)

    t0          = Tint(1);
    tlength     = Tint(2) - Tint(1);
    % NOTE: Does not take leap seconds into account.
    stepSizeSec = nDays*24*60*60;    % seconds.

    dt          = 0:stepSizeSec:tlength;
    out_times   = t0+dt;

end



% Wrapper around solo.db_get_ts() that normalizes the output to a TSeries.
%
% NOTE: solo.db_get_ts() has been observed to return cell array of TSeries
% (instead of TSeries) for Npas, Tpas and Vpas for
% solo.quicklook_main('2020-10-21T00:00:00', '2020-10-28T00:00:00', '/data/solo/data_yuri', ...)
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



% Interpret argument for main function interface. Intended accept and normalize
% arguments which are either
% (1) MATLAB-friendly (numeric/logical), or
% (2) bash script-friendly (strings).
%
function value = interpret_argument_flag(arg)
    assert(isscalar(arg), 'Flag argument is not scalar.')
    
    if isnumeric(arg) || islogical(arg)
        value = logical(arg);
    elseif ischar(arg) && arg=='0'
        value = false;
    elseif ischar(arg) && arg=='1'
        value = true;
    else
        error('Can not interpret argument flag. Illegal format.')
    end
end
