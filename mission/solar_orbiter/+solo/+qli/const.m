%
% Class for collecting constants relating to solo.qli (but not solo.qli.cron
% (almost)).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef const



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % IMPLEMENTATION NOTE: Disabling B (use empty; pretend there is no B data)
    % speeds up solo.qli.quicklooks_24_6_2_h() greatly. Useful for some debugging.
    % Should be enabled by default.
    ENABLE_B = true;

    % Whether to catch plotting exceptions, continue plotting other days/weeks, and
    % then re-raise the last caught exception at the very end. This produces as many
    % quicklooks as possible when one or some quicklooks fail.
    % Should be enabled by default.
    CATCH_PLOT_EXCEPTIONS_ENABLED = true;

    % Path to IRF logo, relative to irfu-matlab root.
    IRF_LOGO_RPATH = 'mission/solar_orbiter/+solo/irf_logo.png';

    % NOTE: Usually found at /data/solo/data_yuri/.
    VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
    VHT_6H_DATA_FILENAME = 'V_RPW.mat';

    % Define on which weekday weekly quicklooks begin
    % -----------------------------------------------
    % NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
    % Therefore using Wednesday as beginning of "week" for weekly quicklooks
    % (until someone complains).
    FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday

    % Host names used for determining whether the processing is "official
    % processing" of quicklooks or not.
    OFFICIAL_PROCESSING_IRFU_HOST_NAMES_CA = {'brain', 'spis'};

    % Directory which shall be used for trying to trigger automounting for
    % official processing.
    OFFICIAL_PROCESSING_AUTOMOUNT_DIR = '/data/solo/';
  end



end
