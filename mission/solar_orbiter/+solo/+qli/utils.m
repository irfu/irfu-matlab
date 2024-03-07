%
% Miscellaneous utility functions. Mostly to collect small shared functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils



  properties(Constant)
    SOAR_URL = 'https://soar.esac.esa.int/';
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Assert that datetime object only contains timestamps which refer to
    % midnight.
    function assert_UTC_midnight_datetime(Dt)
      assert(isa(Dt, 'datetime'))
      assert(strcmp(Dt.TimeZone, 'UTCLeapSeconds'), ...
        'datetime object is not UTC.')
      assert(all(Dt == dateshift(Dt, 'start', 'day'), 'all'), ...
        'datetime object does not only contain timestamps representing midnight.')
    end



    function t = scalar_datetime_to_EpochTT(Dt)
      % NOTE: Might not be the perfect implementation, but it works.

      assert(isa(Dt, 'datetime'))
      assert(strcmp(Dt.TimeZone, 'UTCLeapSeconds'))

      tt2000 = irf.cdf.datevec_to_TT2000(datevec(Dt));
      t = irf.time_array(tt2000);
    end



    % Given an array of datetime representing days to be plotted, derive the
    % corresponding weeks which overlap with the specified days.
    %
    % ARGUMENTS
    % =========
    % DayDtArray
    %       datetime column array. Every timestamp is midnight and represents
    %       the 24h period which begins at that timestamp.
    % firstDayOfWeek
    %       First day of week. datetime convention (1=Sunday, ..., 7=Saturday).
    %
    %
    % RETURN VALUE
    % ============
    % WeekDtArray
    %       datetime column array. Every timestamp is midnight and represents
    %       the week (contiguous 7-day period) which begins at that timestamp.
    %
    function WeekDtArray = derive_weeks(DayDtArray, firstDayOfWeek)
      solo.qli.utils.assert_UTC_midnight_datetime(DayDtArray)
      assert(iscolumn(DayDtArray))

      % Find nearest previous day with specified weekday
      % ------------------------------------------------
      % IMPLEMENTATION NOTE: dateshift(... 'dayofweek' ...) can only search
      % forward. Subtracts days to "round down" instead.
      WeekDtArray = dateshift(DayDtArray - caldays(6), 'dayofweek', firstDayOfWeek);

      % IMPLEMENTATION NOTE: Important to eliminate doubles since, every initial
      % timestamp within the same week will separately generate the same
      % timestamp reresenting the same week.
      WeekDtArray = sort(unique(WeekDtArray));
    end



    % Generate text string with information on data source and when the plot
    % was generated.
    function str = get_data_source_info_string()
      dateStr = char(datetime("now", "Format", "uuuu-MM-dd"));
      str = sprintf( ...
        [ ...
        'Swedish Institute of Space Physics, Uppsala (IRFU), %s.', ...
        ' Data available at %s.' ...
        ], ...
        dateStr, solo.qli.utils.SOAR_URL ...
        );
    end



    % Generate human-readable "context strings" to put at bottom of quicklooks.
    %
    %
    % IMPLEMENTATION NOTES
    % ====================
    % Could be split up into two functions, one per string. It has not been
    % split up both strings can be seen as complementary, or as a top and bottom
    % string.
    % --
    % One can use "R=" etc, but then the text comes too close to the UTC date
    % (reduce font size?).
    %
    %
    % ARGUMENTS
    % =========
    % soloPosTSeries, earthPosTSeries
    %       TSeries with positions for SolO and Earth.
    % Tint
    %       Selected time interval.
    %
    %
    % RETURN VALUES
    % =============
    % Two human-readable one-row (no line feed) strings.
    %       Empty string(s) if no corresponding data for time interval.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function [soloStr, earthStr] = get_context_info_strings(soloPosTSeries, earthPosTSeries, Tint)
      % PROPOSAL: No Tint argument. Caller submits already truncated TSeries.
      %   PRO: One fewer arguments.
      %   CON: Caller has to truncate twice.
      %   CON: Caller might truncate differently for different TSeries.

      assert(isa(soloPosTSeries,  'TSeries'))
      assert(isa(earthPosTSeries, 'TSeries'))
      assert(isa(Tint,            'EpochTT'))
      assert(length(Tint) == 2)

      % NOTE: In principle a lot of execution/time just for obtaining a constant,
      %       but the function is not time-critical so it should not be a problem.
      Units = irf_units;
      AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

      soloPos = soloPosTSeries.tlim(Tint).data;
      if ~isempty(soloPos)
        % NOTE: 2x whitespaces between every value.
        soloStr = sprintf([
          'SolO:', ...
          '  %.2f AU,', ...
          '  EcLat %d\\circ,', ...
          '  EcLon %d\\circ', ...
          ], ...
          soloPos(1,1)/AU_KM, ...
          round(soloPos(1,3)*180/pi), ...
          round(soloPos(1,2)*180/pi) ...
          );
      else
        soloStr = '';
      end

      earthPos = earthPosTSeries.tlim(Tint).data;
      if ~isempty(earthPos)
        earthStr = sprintf('Earth:  EcLon %d\\circ', round(earthPos(1,2)*180/pi));
      else
        earthStr = '';
      end
    end



    % ~Utility function that removes duplicated code from plot functions.
    % NOTE: Function can not simultaneously handle both yyaxis left & right.
    function ensure_axes_data_tick_margins(hAxesArray)
      assert(isa(hAxesArray, 'matlab.graphics.axis.Axes'))

      for i = 1:numel(hAxesArray)
        h = hAxesArray(i);

        h.YLim = solo.qli.ensure_data_tick_margins(...
          h.YTick, [h.YLim(1), h.YLim(2)], h.YScale ...
          );
      end
    end



    function filename = get_plot_filename(Tint)
      assert(isa(Tint, 'EpochTT') && (length(Tint) == 2))

      ett1 = Tint(1);
      utcStr1 = ett1.utc;
      utcStr1 = utcStr1(1:13);
      utcStr1([5,8])=[];

      ett2 = Tint(end);
      utcStr2 = ett2.utc;
      utcStr2 = utcStr2(1:13);
      utcStr2([5,8])=[];

      filename = [utcStr1,'_',utcStr2,'.png'];
    end



    function save_figure_to_file(parentDirPath, Tint)
      % PROPOSAL: Include fig.PaperPositionMode='auto';

      filename = solo.qli.utils.get_plot_filename(Tint);
      filePath = fullfile(parentDirPath, filename);
      print('-dpng', filePath);
    end



    % Simple function for logging number of seconds from previous call.
    % For debugging speed.
    function tBeginSec = log_time(locationStr, tBeginSec)
      tSec = toc(tBeginSec);
      fprintf(1, '%s: %.1f [s]\n', locationStr, tSec)
      tBeginSec = tic();
    end



  end    % methods(Static)



end
