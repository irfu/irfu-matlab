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
    % Whether to enable (permit) having magnetic field data
    % -----------------------------------------------------
    % IMPLEMENTATION NOTE: Disabling B (use empty; pretend there is no B data)
    % speeds up solo.qli.generate_quicklooks_24h_6h_2h() greatly. Useful for
    % some debugging. Should be enabled by default.
    ENABLE_B = true;

    % Whether to enable/disable panels with time-consuming spectra
    % ------------------------------------------------------------
    % Disabling this is useful for debugging and testing (speeds up execution).
    % Should be enabled by default.
    NONWEEKLY_SPECTRA_ENABLED = true;

    % Whether to generate all more than one quicklook of every type (per day)
    % -----------------------------------------------------------------------
    % In practice only affects 6h and 2h quicklooks. Disabling this is useful
    % for debugging and testing (speeds up execution). Should be enabled by
    % default.
    NONWEEKLY_ALL_PLOTS_ENABLED = true;
    %NONWEEKLY_ALL_PLOTS_ENABLED = false;



    % Whether to catch plotting exceptions, continue plotting other days/weeks,
    % and then re-raise the last caught exception at the very end. This produces
    % as many quicklooks as possible when one or some quicklooks fail. Should be
    % enabled by default.
    CATCH_PLOT_EXCEPTIONS_ENABLED = true;

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