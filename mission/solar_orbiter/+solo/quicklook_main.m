%
% Generate multiple IRFU-local types of quicklooks (files) for SolO data (not
% just RPW) for a specified time interval.
%
% Intended for batch processing.
%
%
% ARGUMENTS
% =========
% utcBegin, utcEnd : Strings.
%       Defines time interval for which quicklooks should be generated.
% outputDir : Path to directory.
%       NOTE: Will create subdirectories if not pre-existing.
%
%
% NOTE: Uses SPICE implicitly, and therefore relies on some path convention.
% NOTE: Still uses some hard-coded paths for locating datasets.
%
%
% Initially created ~<2021-03-11, based on code by Konrad Steinvall, IRF,
% Uppsala, Sweden. Modified by Erik P G Johansson.
%
function quicklook_main(utcBegin, utcEnd, speedDataDir, outputDir)

% IMPLEMENTATION NOTE: Needed to make "DB" work.
irf



% Specify Solar Orbiter database.
% For data on the server the following two are sufficient. You also need
% V_RPW.mat and V_RPW_1h.mat, found on solo/data_yuri.
solo.db_init('local_file_db', '/data/solo/');
solo.db_init('local_file_db', '/data/solo/data_irfu');

% Setup cache
solo.db_init('db_cache_size_max',4096)
solo.db_cache('on','save')

% Speeds up solo.quicklooks_24_6_2_h(). Useful for testing.
DISABLE_B = 0;

% Specify subdirectories for saving the respective types of plots.
PATHS.path_2h  = fullfile(outputDir, '2h' );
PATHS.path_6h  = fullfile(outputDir, '6h' );
PATHS.path_24h = fullfile(outputDir, '24h');
PATHS.path_1w  = fullfile(outputDir, '1w' );



% Specify time interval for plots
TimeInterval = irf.tint(utcBegin, utcEnd);

% Create subdirectories
% NOTE: Works if subdirectories pre-exist.
for fnCa = fieldnames(PATHS)'
    dirPath = PATHS.(fnCa{1});
    [parentDir, dirBasename, dirSuffix] = fileparts(dirPath);
    [success, msg] = mkdir(parentDir, [dirBasename, dirSuffix]);
    assert(success, 'Failed to create directory "%s". %s', dirPath, msg)
end



%=============================================
% Run the code for 2-, 6-, 24-hour quicklooks
%=============================================
times_1d = make_tints(TimeInterval,1);% Daily time-intervals

% Load data
% This is the .mat file containing RPW speeds at 1h resolution.
% The file should be in the current path. This file can be found in
% /solo/data/data_yuri/.
S = load(fullfile(speedDataDir, 'V_RPW_1h'));

for iTint=1:length(times_1d)-1
    % Time interval
    Tint=irf.tint(times_1d(iTint),times_1d(iTint+1));

    Data.Vrpw = S.V_RPW_1h.tlim(Tint);

    % E-field:
    Data.E = db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    Data.Ne = db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    Data.B = db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);

    % Proton & alpha temperature:
    Data.Tpas = db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    Data.Vpas = db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);

    % Proton & alpha density:
    Data.Npas = db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

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



%===================================
% Run the code for weekly overviews
%===================================
times_7d = make_tints(TimeInterval,7);% weekly time-intervals


% Load data
% This is the .mat file containing RPW speeds at 6h resolution.
% The file should be in the same folder as this script (quicklook_main).
S = load(fullfile(speedDataDir, 'V_RPW'));

for iTint=1:length(times_7d)-1
    % Time interval
    Tint = irf.tint(times_7d(iTint),times_7d(iTint+1));
    data2.Vrpw = S.V_RPW.tlim(Tint);

    % E-field:
    data2.E = db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);

    % RPW density:
    data2.Ne = db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);

    % B-field:
    data2.B = db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);

    % Proton & alpha temperature:
    data2.Tpas = db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);

    % Proton & alpha velocity:
    data2.Vpas = db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);

    % Proton & alpha density:
    data2.Npas = db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);

    % Solar Orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".
    data2.solopos = solo.get_position(Tint,'frame','SOLO_SUN_RTN');

    % Plot data and save figure
    solo.quicklooks_7days(data2,PATHS,Tint)

end



% Auxilliary function
%
function out_times=make_tints(inputTint,days)

    t0       = inputTint(1);
    tlength  = inputTint(2)-inputTint(1);
    stepsize = days*24*60*60; %seconds

    dt = 0:stepsize:tlength;
    out_times = t0+dt;

end



end    % function

% Wrapper around solo.db_get_ts() that normalizes the output to a TSeries.
%
% NOTE: solo.db_get_ts() has been observed to return cell array of TSeries
% (instead of TSeries) for Npas, Tpas and Vpas for
% solo.quicklook_main('2020-10-21T00:00:00', '2020-10-28T00:00:00', '/data/solo/data_yuri', ...)
% There might be other cases but those are as of yet unknown.
% /Erik P G Johansson 2021-03-22
%
function Ts = db_get_ts(varargin)
    Ts = solo.db_get_ts(varargin{:});
    
    % Normalize (TSeries or cell array) --> TSeries.
    if iscell(Ts)
        Ts = cell_array_TS_to_TS(Ts);
    end
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
