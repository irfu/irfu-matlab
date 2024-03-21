%
% Class for collecting constants relating to solo.qli (but not solo.qli.offgen
% (almost)).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef const
  % PROPOSAL: Move OFFICIAL_GENERATION_IRFU_HOST_NAMES_CA and
  %           OFFICIAL_GENERATION_AUTOMOUNT_DIR to solo.qli.offgen somehow.



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
    ENABLE_B = true;    % DEFAULT
    %ENABLE_B = false;

    % Whether to enable/disable panels with time-consuming spectra
    % ------------------------------------------------------------
    % Disabling this is useful for debugging and testing (speeds up execution).
    % Should be enabled by default.
    NONWEEKLY_SPECTRA_ENABLED = true;    % DEFAULT

    % Whether to generate more than one quicklook (file) of every type (per day)
    % --------------------------------------------------------------------------
    % Only affects 6h and 2h quicklooks in practice. Disabling this is useful
    % for debugging and testing (speeds up execution). Should be enabled by
    % default.
    NONWEEKLY_ALL_PLOTS_ENABLED = true;    % DEFAULT
    %NONWEEKLY_ALL_PLOTS_ENABLED = false;



    % Whether to catch plotting exceptions, continue plotting other days/weeks,
    % and then re-raise the last caught exception at the very end. This produces
    % as many quicklooks as possible when one or some quicklooks fail. Should be
    % enabled by default.
    CATCH_PLOT_EXCEPTIONS_ENABLED = true;

    % NOTE: Below files are normally found at brain:/data/solo/data_yuri/.
    VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
    VHT_6H_DATA_FILENAME = 'V_RPW.mat';

    % Define on which weekday weekly quicklooks begin
    % -----------------------------------------------
    % NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
    % Therefore using Wednesday as beginning of "week" for weekly quicklooks
    % (until someone complains).
    FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday

    % Host names used for determining whether the current generation of
    % quicklooks is "official generation" of quicklooks or not.
    OFFICIAL_GENERATION_IRFU_HOST_NAMES_CA = {'brain', 'spis'};

    % Directory which shall be used for trying to trigger automounting for
    % official generation of quicklooks.
    OFFICIAL_GENERATION_AUTOMOUNT_DIR = '/data/solo/';
  end



end
