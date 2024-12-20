%
% Class for collecting constants relating to solo.qli.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef const
  % PROPOSAL: Move VHT_1H_DATA_FILENAME, VHT_6H_DATA_FILENAME to solo.qli.batch
  %           somehow.
  %   PROBLEM: generate_quicklooks_*_using_DB_SPICE() use them.
  %     PROPOSAL: Also hardcode values in demo file.
  %
  % PROPOSAL: Replace B_SPECTRA_ENABLED, NONWEEKLY_ALL_PLOTS_ENABLED,
  %           NONWEEKLY_ALL_PLOTS_ENABLED
  %           with optional arguments, "settings" (varargin).
  %   NOTE: Would need to expose the same settings arguments in
  %         solo.qli.batch.generate_quicklooks().
  %   PRO: Calling MTEST code could set values (override defaults).
  %
  % PROPOSAL: Abolish CATCH_PLOT_EXCEPTIONS_ENABLED.
  %
  % PROBLEM: "EMPTY_DT_ARRAY" does not imply that it is UTCLeapSeconds.
  %   PROPOSAL: Abbreviation for UTC datetime objects/arrays.
  %     NOTE: Already exists UMD=UTC Midnight Date. Not defined to be datetime
  %           specifically, though seems to be datetime in practice.
  %     PROPOSAL: UDT  = UTC datetime
  %     PROPOSAL: UMDT = UTC midnight datetime
  %       PRO: Fits better with simultaneously using MDT = Midnight datetime.
  %     PROPOSAL: MUDT = Midnight UTC datetime
  %       PRO: Fits better with simultaneously using UDT.
  %     --
  %     CON: There is value in separately specifying (1) the properties of
  %          timestamps, (2) the time format of the timestamp, and (3) what a
  %          timestamp signifies.
  %       Ex:
  %         Properties:   Midnight, UTC
  %         Time formats: datetime, SDN, date vector
  %         Other:        file modification date, ~dataset date
  %       CON: In practice, only datetime is used for QLI, and should be
  %            converged upon for other applications.
  %     CON: More natural to describe "algorithms" in terms of abbreviations
  %          which do not specify time format, e.g. FMD.
  %     PRO: Shortens variables.
  %       Ex: UmdDt-->Umdt, FmdDt-->Fmdt



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant)

    %========================================
    % Enable/disable functionality for speed
    %========================================

    % Whether to enable two time-consuming spectra based on B (24h6h2h)
    % -----------------------------------------------------------------
    % IMPLEMENTATION NOTE: Disabling these speeds up
    % solo.qli.generate_quicklooks_24h_6h_2h() greatly. Useful for some
    % debugging. Should be enabled by default.
    B_SPECTRA_ENABLED = true;    % DEFAULT
    %B_SPECTRA_ENABLED = false;

    % Whether to generate any 6h and 2h quicklooks (files)
    % ----------------------------------------------------
    % Disabling this is useful for debugging and testing (speeds up execution).
    % Should be enabled by default.
    NONWEEKLY_6H_2H_PLOTS_ENABLED = true;
    %NONWEEKLY_6H_2H_PLOTS_ENABLED = false;

    % Whether to generate more than one quicklook (file) of every type (per day)
    % --------------------------------------------------------------------------
    % Only affects 6h and 2h quicklooks in practice. Disabling this is useful
    % for debugging and testing (speeds up execution). Should be enabled by
    % default.
    NONWEEKLY_ALL_PLOTS_ENABLED = true;    % DEFAULT
    %NONWEEKLY_ALL_PLOTS_ENABLED = false;



    %===============
    % Miscellaneous
    %===============

    % Whether to catch plotting exceptions, continue plotting other days/weeks,
    % and then re-raise the last caught exception at the very end. This produces
    % as many quicklooks as possible when one or some quicklooks fail. Should be
    % enabled by default.
    CATCH_PLOT_EXCEPTIONS_ENABLED = true;

    % Log level used for irf.log() when logging exceptions which are caught due
    % to CATCH_PLOT_EXCEPTIONS_ENABLED=true.
    LOG_LEVEL_CAUGHT_EXCEPTIONS = 'C';

    % NOTE: Below files are normally found at brain:/data/solo/data_yuri/.
    VHT_1H_DATA_FILENAME = 'V_RPW_1h.mat';
    VHT_6H_DATA_FILENAME = 'V_RPW.mat';

    % Define on which weekday weekly quicklooks begin
    % -----------------------------------------------
    % NOTE: First day of data (launch+2 days) is 2020-02-12, a Wednesday.
    % Therefore using Wednesday as beginning of "week" for weekly quicklooks
    % (until someone complains).
    FIRST_DAY_OF_WEEK = 4;   % 2 = Monday; 4 = Wednesday

    % Constant which is useful to have for defining tests.
    EMPTY_DT_ARRAY = NaT([0, 1], 'TimeZone', 'UTCLeapSeconds');

    % NOTE: Used for putting in human-readable strings (on QLIs), not for
    % loading data over the internet.
    SOAR_URL = 'https://soar.esac.esa.int/';

  end



end
