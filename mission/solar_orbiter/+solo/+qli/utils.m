%
% Miscellaneous utility functions. Mostly to collect small shared functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Given an array of datetime representing days to be plotted, derive the
    % corresponding weeks which overlap with the specified days.
    %
    % ARGUMENTS
    % =========
    % UmdDtArray
    %       UMD datetime column array. Every timestamp represents
    %       the 24h period which begins at that timestamp.
    % firstDayOfWeek
    %       First day of week. datetime convention (1=Sunday, ..., 7=Saturday).
    %
    %
    % RETURN VALUE
    % ============
    % WeekUmdDtArray
    %       UMD datetime column array. Every timestamp represents the week
    %       (contiguous 7-day period) which begins at that timestamp.
    %
    function WeekUmdDtArray = derive_weeks(UmdDtArray, firstDayOfWeek)
      irf.dt.assert_UTC_midnight(UmdDtArray)
      assert(iscolumn(UmdDtArray))

      % Find nearest previous day with specified weekday
      % ------------------------------------------------
      % IMPLEMENTATION NOTE: dateshift(... 'dayofweek' ...) can only search
      % forward. Subtracts days to "round down" instead.
      WeekUmdDtArray = dateshift(UmdDtArray - caldays(6), 'dayofweek', firstDayOfWeek);

      % IMPLEMENTATION NOTE: Important to eliminate doubles since, every initial
      % timestamp within the same week will separately generate the same
      % timestamp reresenting the same week.
      WeekUmdDtArray = sort(unique(WeekUmdDtArray));
    end



    % Specify and create subdirectories to place quicklooks in.
    function OutputPaths = create_output_directories(outputDir)
      % PROBLEM: Hardcoded subdirectory names in code.
      % IMPLEMENTATION NOTE: Function is useful partly since it can be shared
      % with test code.

      OutputPaths.dir2h  = fullfile(outputDir, '2h' );
      OutputPaths.dir6h  = fullfile(outputDir, '6h' );
      OutputPaths.dir24h = fullfile(outputDir, '24h');
      OutputPaths.dir1w  = fullfile(outputDir, '1w' );

      %=========================================================
      % Create subdirectories for different types of quicklooks
      %=========================================================
      % NOTE: Works (without warnings) also if subdirectories pre-exist
      % NOTE: Creates subdirectories for all four quicklook types, regardless of
      %       whether they will actually be generated or not.
      %   NOTE: This simplifies bash wrapper scripts that copy content of
      %         sub-directories to analogous sub-directories, also when not all
      %         quicklook types are generated.
      for fnCa = fieldnames(OutputPaths)'
        dirPath = OutputPaths.(fnCa{1});
        if ~exist(dirPath, 'dir')
          [success, msg] = mkdir(dirPath);
          assert(success, 'Failed to create directory "%s": "%s"', dirPath, msg)
        end
      end
    end



    % Read ONE zVariable for *constant* metadata which can not be loaded using
    % solo.db_get_ts() due to being constant CDF metadata(?), e.g. not having
    % DEPEND_0..
    %
    %
    % NOTE: Will only load the first dataset found since the zVariable is
    %       assumed to be constant across datasets.
    %
    %
    % EXAMPLES OF zVARIABLES WHICH CAN BE LOADED WITH THIS FUNCTION, BUT NOT
    % WITH solo.db_get_ts()
    % ======================================================================
    % Ex: solo_L2_swa-pas-eflux + zVariable "Energy"
    %     solo.db_get_ts() yields error message: "Data does not contain DEPEND_0
    %     or DATA"
    %
    % Ex: solo_L2_rpw-tnr-surv-cdag + zVariable "TNR_BAND_FREQ"
    %     solo.get_db_ts() function seems to fail to create the TSeries object
    %     because the DEPEND_0 zVariable is of different size compared to the
    %     requested zVariable.
    %     solo.db_get_ts() yields error message:
    %     """"
    %     Output argument "TS" (and maybe others) not assigned during call to "irf.ts_scalar".
    %     Error in solo.variable2ts (line 44)
    %     ts = feval(['irf.ts_' varType],v.DEPEND_0.data,data);
    %     Error in solo_db/get_ts (line 327)
    %             res = solo.variable2ts(v);
    %     Error in solo.db_get_ts (line 50)
    %     res = SOLO_DB.get_ts(filePrefix,varName,tint);
    %     """"
    %
    %
    % RETURN VALUE
    % ============
    % zvData
    %       Array, if found at least one dataset.
    %       [], if no matching dataset was found.
    %
    function zvData = read_constant_metadata(filePrefix, zvName, Tint)
      % PROPOSAL: Refactor to use dataobj().

      FileArray = solo.db_list_files(filePrefix, Tint);
      if ~isempty(FileArray)
        FileArray(1);
        filePath = fullfile(FileArray(1).path, FileArray(1).name);

        % NOTE: Reads CDFs using cdfread() which is a MATLAB function (i.e. not
        %       dataobj(), not NASA SPDF). *Might* be faster (or might not)
        %       since only reading one specific zVariable.
        zvCa   = cdfread(filePath, 'variables', zvName);
        zvData = zvCa{1};
      else
        zvData = [];
      end
    end



    % Wrapper around solo.db_get_ts() which normalizes the output to always
    % return one TSeries object, or [].
    %
    % NOTE: solo.db_get_ts() returns a cell array of TSeries instead of a single
    % TSeries when the underlying code thinks that the underlying CDFs do not
    % have consistent metadata. See solo.db_get_ts().
    %
    function Ts = db_get_ts(varargin)
      % TODO-NI: Document example for which solo.db_get_ts() returns cell array.

      temp = solo.db_get_ts(varargin{:});

      % Normalize (TSeries or cell array) --> TSeries.
      if iscell(temp)
        TsCa   = temp;   % Rename variable.
        nCells = numel(TsCa);
        Ts     = TsCa{1};

        if nCells>1
          for iCell = 2:nCells    % NOTE: Begins at 2.
            Ts = Ts.combine(TsCa{iCell});
          end
        end
      else
        Ts = temp;
      end

    end



    % Use SPICE to get Solar Orbiter's position as
    % [soloSunDistance, soloEclLongitude, soloEclLatitude].
    % Longitude and latitude are in radians.
    %
    % NOTE: Uses SPICE kernels and "solo.get_position()" indirectly through
    % irfu-matlab which can itself load SPICE kernels(!).
    function soloPosRadLonLatTSeries = get_SolO_position(Tint)
      assert((length(Tint) == 2) & isa(Tint, 'EpochTT'))

      % See solo.qli.utils.get_Earth_position() (in this file) for information
      % on the coordinate system.
      % NOTE: Function automatically returns data with a sampling rate of
      % 1 data point/hour.
      soloPosXyz = solo.get_position(Tint, 'frame', 'ECLIPJ2000');

      if ~isempty(soloPosXyz)
        [soloSunDistance, soloEclLongitudeRad, soloEclLatitudeRad] = cspice_reclat(soloPosXyz.data');
        soloPosRadLonLatTSeries = irf.ts_vec_xyz(soloPosXyz.time, ...
          [soloSunDistance', soloEclLongitudeRad', soloEclLatitudeRad']);
      else
        %soloPosRadLonLat = soloPosXyz;
        soloPosRadLonLatTSeries = TSeries();   % Empty TSeries.
      end
    end



    % Use SPICE to get Earth's position as
    % [earthSunDistance, earthEclLongitude, earthEclLatitude].
    % Longitude and latitude are in radians.
    %
    function earthPosRadLonLatTSeries = get_Earth_position(Tint, dtSec)
      %=========================================================================
      % Arguments for cspice_spkpos()
      % -----------------------------
      % 17  ECLIPJ2000  Ecliptic coordinates based upon the
      %                 J2000 frame.
      %
      %                 The value for the obliquity of the
      %                 ecliptic at J2000 is taken from page 114
      %                 of [7] equation 3.222-1. This agrees with the
      %                 expression given in [5].
      %
      % Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
      % --
      % 'LT+S'     Correct for one-way light time and
      %            stellar aberration using a Newtonian
      %            formulation. This option modifies the
      %            position obtained with the 'LT' option
      %            to account for the observer's velocity
      %            relative to the solar system
      %            barycenter. The result is the apparent
      %            position of the target---the position
      %            as seen by the observer.
      %
      % Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_spkpos.html
      %=========================================================================
      assert((length(Tint) == 2) & isa(Tint, 'EpochTT'))
      assert(isnumeric(dtSec))

      et = Tint.start.tts : dtSec : Tint.stop.tts;

      earthPosXyz = cspice_spkpos('Earth', et, 'ECLIPJ2000', 'LT+s', 'Sun');

      if ~isempty(earthPosXyz)
        [earthSunDistance, earthEclLongitudeRad, earthEclLatitudeRad] = cspice_reclat(earthPosXyz);
        earthPos = [earthSunDistance', earthEclLongitudeRad', earthEclLatitudeRad'];

        Tlength = Tint(end)-Tint(1);
        dTimes  = 0:dtSec:Tlength;
        Times   = Tint(1)+dTimes;
        earthPosRadLonLatTSeries = irf.ts_vec_xyz(Times, earthPos);
      else
        earthPosRadLonLatTSeries = TSeries();   % Empty TSeries.
      end
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
        dateStr, solo.qli.const.SOAR_URL ...
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
    function [soloStr, earthStr] = get_context_info_strings(...
        soloPosTSeries, earthPosTSeries, Tint)
      % PROPOSAL: No Tint argument. Caller submits already truncated TSeries.
      %   PRO: One fewer arguments.
      %   CON: Caller has to truncate twice.
      %   CON: Caller might truncate differently for different TSeries.

      assert(isa(soloPosTSeries,  'TSeries'))
      assert(isa(earthPosTSeries, 'TSeries'))
      assert(isa(Tint,            'EpochTT'))
      assert(length(Tint) == 2)

      %======
      % SolO
      %======
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

      %=======
      % Earth
      %=======
      earthPos = earthPosTSeries.tlim(Tint).data;
      if ~isempty(earthPos)
        earthStr = sprintf('Earth:  EcLon %d\\circ', round(earthPos(1,2)*180/pi));
      else
        earthStr = '';
      end
    end



    % ~Utility function which removes duplicated code from plot functions.
    %
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



    % For debugging.
    function print_Y_axis_info(h)
      % PROPOSAL: Use irf.log('debug', ...)

      fprintf('h.YLimMode  = "%s"\n', h.YLimMode)
      fprintf('h.YTickMode = "%s"\n', h.YLimMode)
      fprintf('h.YLim  = [%s]\n', sprintf('%f  ', h.YLim))
      fprintf('h.YTick = [%s]\n', sprintf('%f  ', h.YTick))
    end



    % (1) Set YLim automatically (from data; with margins)
    % (2) Set YTick automatically.
    %
    % NOTE: Function can not simultaneously handle both yyaxis left & right.
    % NOTE: MATLAB's automatic setting of y ticks for log scale (and which is used)
    %       can be bad. May therefore want to set y ticks for panels with log scale.
    function set_YLim_YTick_automatically(hAxesArray)
      assert(isa(hAxesArray, 'matlab.graphics.axis.Axes'))

      %h = hAxesArray(1);
      %solo.qli.utils.print_Y_axis_info(h)

      %=======================================================================
      % Automatically set preliminary YLim (y limits) and final YTick (y tick
      % positions) for selected axes.
      %=======================================================================
      % Set axes y range (YLim) to only cover the data (plus rounding outwards to
      % ticks).
      set(hAxesArray, 'YLimMode', 'auto')
      % Auto-generate ticks (YTick; y values at which there should be ticks).
      set(hAxesArray, 'YTickMode', 'auto')

      %---------------------------------------------------------------------------
      % IMPORTANT: Read YLim without using the return result ("do nothing")
      % -------------------------------------------------------------------
      % IMPLEMENTATION NOTE: THIS COMMAND SHOULD THEORETICALLY NOT BE NEEDED,
      % BUT IS NEEDED FOR THE YLim VALUES TO BE SET PROPERLY. MATLAB BUG?!
      % This behaviour has been observed on Erik P G Johansson's laptop
      % "irony" (MATLAB R2019b, Ubuntu Linux) as of 2023-05-25.
      % Ex: (Re-)scaling of panel 5, 2022-02-23T10-12 (2h plot).
      % 2024-03-22: Refactored the code. Not sure if needed.
      get(hAxesArray, 'YLim');
      %---------------------------------------------------------------------------
      % Prevent later setting of YLim (next command) from generating new ticks.
      set(hAxesArray, 'YTickMode', 'manual')

      %solo.qli.utils.print_Y_axis_info(h)
      solo.qli.utils.ensure_axes_data_tick_margins(hAxesArray)
      %solo.qli.utils.print_Y_axis_info(h)
    end



    % (1) Set YLim automatically (from data; with margins),
    % (2) Keep YTick as is.
    %
    % NOTE: See set_YLim_YTick_automatically().
    function set_YLim_automatically(hAxesArray)
      assert(isa(hAxesArray, 'matlab.graphics.axis.Axes'))

      %h = hAxesArray(1);
      %solo.qli.utils.print_Y_axis_info(h)

      %=========================================================================
      % Automatically set YLim (y limits) but keep old YTick (y tick positions)
      % for selected axes.
      %=========================================================================
      set(hAxesArray, 'YTickMode', 'manual')
      %get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?
      set(hAxesArray, 'YLimMode',  'auto')
      %get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?

      %solo.qli.utils.print_Y_axis_info(h)
      solo.qli.utils.ensure_axes_data_tick_margins(hAxesArray)
      %solo.qli.utils.print_Y_axis_info(h)
    end



    % Ex: 20240313T20_20240313T22.png
    function filename = create_quicklook_filename(Tint)
      assert(isa(Tint, 'EpochTT') && (length(Tint) == 2))

      ett1           = Tint(1);
      utcStr1        = ett1.utc;
      utcStr1        = utcStr1(1:13);
      utcStr1([5,8]) = [];

      ett2           = Tint(2);
      utcStr2        = ett2.utc;
      utcStr2        = utcStr2(1:13);
      utcStr2([5,8]) = [];

      filename = [utcStr1, '_', utcStr2,'.png'];
    end



    % Parse a quicklook filename.
    %
    % ARGUMENTS
    % =========
    % filename
    %       Quicklook filename created with
    %       solo.qli.utils.create_quicklook_filename().
    %
    % RETURN VALUES
    % =============
    % Dt1, Dt2
    %       datetime. Beginning and end timestamps, if the filename is a valid
    %       quicklook filename.
    %       [], if non-quicklook filename.
    function [Dt1, Dt2] = parse_quicklook_filename(filename)
      TIMESTAMP_RE_CA = {'20[0-9]{2}', '[01][0-9]', '[0-3][0-9]', 'T', '[0-2][0-9]'};

      % Ex: 20240313T20_20240313T22.png
      [subStrCa, ~, isPerfectMatch] = irf.str.regexp_str_parts(...
        filename, {TIMESTAMP_RE_CA{:}, '_', TIMESTAMP_RE_CA{:}, '.png'}, 'PERMIT_NON_MATCH');

      if isPerfectMatch
        year1  = str2double(subStrCa{1});
        month1 = str2double(subStrCa{2});
        day1   = str2double(subStrCa{3});
        hour1  = str2double(subStrCa{4+1});

        year2  = str2double(subStrCa{6+1});
        month2 = str2double(subStrCa{6+2});
        day2   = str2double(subStrCa{6+3});
        hour2  = str2double(subStrCa{6+4+1});

        DtArray = datetime(...
          [year1, month1, day1, hour1, 0, 0;
          year2, month2, day2, hour2, 0, 0], ...
          'TimeZone', 'UTCLeapSeconds');
        Dt1 = DtArray(1);
        Dt2 = DtArray(2);
      else
        Dt1 = [];
        Dt2 = [];
      end
    end



    function save_figure_to_file(parentDirPath, Tint)
      % PROPOSAL: Include fig.PaperPositionMode='auto';

      filename = solo.qli.utils.create_quicklook_filename(Tint);
      filePath = fullfile(parentDirPath, filename);
      print('-dpng', filePath);
    end



    % Simple function for logging number of seconds from previous call.
    % For debugging speed of execution.
    function tBeginSec = log_time(locationStr, tBeginSec)
      tSec = toc(tBeginSec);
      %fprintf(1, '%s: %.1f [s]\n', locationStr, tSec)
      irf.log('n', sprintf('%s: %.1f [s]', locationStr, tSec))
      tBeginSec = tic();
    end



  end    % methods(Static)



end
