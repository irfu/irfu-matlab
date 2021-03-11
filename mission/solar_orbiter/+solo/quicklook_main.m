%% Specify solar orbiter database.
% For data on the server the following two are sufficient. You also need
% V_RPW.mat and V_RPW_1h.mat, found on solo/data_yuri.

solo.db_init('local_file_db','/Volumes/solo/');
solo.db_init('local_file_db','/Volumes/solo/data_irfu');

% Setup cache
solo.db_init('db_cache_size_max',4096)
solo.db_cache('on','save')
%% Specify time interval for plots
plot_interval = irf.tint('2020-06-01T00:00:00.00Z','2020-08-01T00:00:00.00Z');


%% Specify folders for saving the plots

paths.path_2h=[pwd,'/2h']; %Path to folder for 2-hour overviews
paths.path_6h=[pwd,'/6h']; %Path to folder for 6-hour overviews
paths.path_24h=[pwd,'/24h']; %Path to folder for 24-hour overviews
paths.path_1w=[pwd,'/1w']; %Path to folder for 1w overviews


%% Run the code for 2, 6, 24 hours.
times_1d=make_tints(plot_interval,1);% Daily time-intervals

for iTint=1:length(times_1d)-1
    %Time interval
    Tint=irf.tint(times_1d(iTint),times_1d(iTint+1));
    
    %Load data
    load('V_RPW_1h'); %This is the .mat file containing RPW speeds at 1h resolution.
    % The file should be in the current path. This file can be found in
    % solo/data_yuri
    
    data.Vrpw = V_RPW_1h.tlim(Tint);
    
    %E-field:
    data.E = solo.db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);
    
    %RPW denisty:
    data.Ne = solo.db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);
    
    %B-field:
    data.B = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);

    %Proton & alpha temperature:
    data.Tpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);   
    
    %Proton & alpha velocity:
    data.Vpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);
    
    %Proton & alpha density:
    data.Npas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);    
    
    %Solar orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".    
    data.solopos = solo.get_position(Tint,'frame','SOLO_SUN_RTN');       

    %Plot data and save figure
    solo.quicklooks_24_6_2_h(data,paths,Tint)
     
end


%% Run the code for weekly overviews.
times_7d=make_tints(plot_interval,7);% weekly time-intervals

for iTint=1:length(times_7d)-1
    %Time interval
    Tint=irf.tint(times_7d(iTint),times_7d(iTint+1));
    
    %Load data
    load('V_RPW'); %This is the .mat file containing RPW speeds at 6h resolution.
    % The file should be in the same folder as this script (quicklook_main)
    data2.Vrpw = V_RPW.tlim(Tint);
    
    %E-field:
    data2.E = solo.db_get_ts('solo_L3_rpw-bia-efield-10-seconds', 'EDC_SRF', Tint);
    
    %RPW denisty:
    data2.Ne = solo.db_get_ts('solo_L3_rpw-bia-density-10-seconds', 'DENSITY', Tint);
    
    %B-field:
    data2.B = solo.db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);
    
    %Proton & alpha temperature:
    data2.Tpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','T', Tint);   
    
    %Proton & alpha velocity:
    data2.Vpas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','V_SRF', Tint);
    
    %Proton & alpha density:
    data2.Npas = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N', Tint);    
    
    %Solar orbiter position
    % Note: solopos uses SPICE, but should be taken care of by
    % "solo.get_position".    
    data2.solopos = solo.get_position(Tint,'frame','SOLO_SUN_RTN');       

    %Plot data and save figure
    solo.quicklooks_7days(data2,paths,Tint)
     
end



%% Auxilliary function

function out_times=make_tints(inputTint,days)

t0=inputTint(1);
tlength = inputTint(2)-inputTint(1);
stepsize = days*24*60*60; %seconds

dt = 0:stepsize:tlength;
out_times = t0+dt;

end

