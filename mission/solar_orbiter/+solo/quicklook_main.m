%% Specify Solar Orbiter database.
% For data on the server the following two are sufficient. You also need
% V_RPW.mat and V_RPW_1h.mat, found on solo/data_yuri.

%clear all   % TEMP

solo.db_init('local_file_db','/data/solo/');
solo.db_init('local_file_db','/data/solo/data_irfu');

% Setup cache
solo.db_init('db_cache_size_max',4096)
solo.db_cache('on','save')



%% Specify time interval for plots
TIME_INTERVAL = irf.tint('2020-06-01T00:00:00.00Z','2020-06-02T00:00:00.00Z');

%OUTPUT_DIR = '/home/erjo/temp/qli';
OUTPUT_DIR = '.';
SPEED_DIR  = '/data/solo/data_yuri';

% Speeds up solo.quicklooks_24_6_2_h which is useful for testing purposes.
DISABLE_B = 0;

%% Specify folders for saving the plots

PATHS.path_2h  = fullfile(OUTPUT_DIR, '2h' );    % Path to folder for 2-hour overviews
PATHS.path_6h  = fullfile(OUTPUT_DIR, '6h' );    % Path to folder for 6-hour overviews
PATHS.path_24h = fullfile(OUTPUT_DIR, '24h');    % Path to folder for 24-hour overviews
PATHS.path_1w  = fullfile(OUTPUT_DIR, '1w' );    % Path to folder for 1w overviews



%% Run the code for 2, 6, 24 hours.
times_1d = make_tints(TIME_INTERVAL,1);% Daily time-intervals

% Load data
% This is the .mat file containing RPW speeds at 1h resolution.
% The file should be in the current path. This file can be found in
% /solo/data/data_yuri/.
S = load(fullfile(SPEED_DIR, 'V_RPW_1h'));

for iTint=1:length(times_1d)-1
    % Time interval
    Tint=irf.tint(times_1d(iTint),times_1d(iTint+1));

    Data.Vrpw = S.V_RPW_1h.tlim(Tint);

    % E-field:
    Data.E = solo.db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    Data.Ne = solo.db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    Data.B = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);

    % Proton & alpha temperature:
    Data.Tpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    Data.Vpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);

    % Proton & alpha density:
    Data.Npas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

    % Solar Orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".
    Data.solopos = solo.get_position(Tint,'frame','SOLO_SUN_RTN');

    if DISABLE_B
        Data.B = [];
    end
    %Plot data and save figure
    solo.quicklooks_24_6_2_h(Data,PATHS,Tint)

end

%return    % TEST

%% Run the code for weekly overviews.
times_7d = make_tints(TIME_INTERVAL,7);% weekly time-intervals


% Load data
% This is the .mat file containing RPW speeds at 6h resolution.
% The file should be in the same folder as this script (quicklook_main).
S = load(fullfile(SPEED_DIR, 'V_RPW'));

for iTint=1:length(times_7d)-1
    % Time interval
    Tint = irf.tint(times_7d(iTint),times_7d(iTint+1));
    data2.Vrpw = S.V_RPW.tlim(Tint);

    % E-field:
    data2.E = solo.db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    data2.Ne = solo.db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    data2.B = solo.db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);

    % Proton & alpha temperature:
    data2.Tpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    data2.Vpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);

    % Proton & alpha density:
    data2.Npas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

    % Solar Orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".
    data2.solopos = solo.get_position(Tint,'frame','SOLO_SUN_RTN');

    % Plot data and save figure
    solo.quicklooks_7days(data2,PATHS,Tint)

end



%% Auxilliary function

function out_times=make_tints(inputTint,days)

    t0       = inputTint(1);
    tlength  = inputTint(2)-inputTint(1);
    stepsize = days*24*60*60; %seconds

    dt = 0:stepsize:tlength;
    out_times = t0+dt;

end
