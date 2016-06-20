% SCRIPT_ORBIT_COVERAGE.m
%
% % 1: No need to run. 
% compare_different_years_of_omni_data
%
% % 2: Define regions of interest and see how much time is spent within
%      each bin. Plots regions in 2D-plane.
% region_coverage
%
% % 3: Plots axisymmetric 3D regions of interest.
% regions_of_interest

%% compare_different_years_of_omni_data.m
% Load omni data for 3 different years and plot some reference bowshock and
% magnetopause locations.
compare_different_years_of_omni_data


%% region_coverage.m
% Load orbit data, loads omni data to help define regions of interest, bin
% the time spent in each region.
region_coverage

%% regions_of_interest.m
% Plots axisymmetric regions of interest.
useMemoValues = 0; % uses values from region_coverage - script
regions_of_interest