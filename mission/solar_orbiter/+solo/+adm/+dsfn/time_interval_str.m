%
% Collections of (two) functions for creating and parsing time interval
% strings for supporting solo.adm.dsfn.DatasetFilename.
%
% NOTE: Tested directly by solo.adm.dsfn.time_interval_str___UTEST and
% indirectly by solo.adm.dsfn.DatasetFilename___UTEST.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef time_interval_str
% PROPOSAL: Return special value(s) for parsing illegal string.
%   PRO: Useful for rigorously distinguishing datasets/non-datasets.
%   PROPOSAL: timeIntervalFormat == special value, e.g. [].
%   PROPOSAL: Special return value flag.
%
% PROPOSAL: New abbrev. TIFID = Time Interval Format ID
%
% PROPOSAL: Class name?
%     tis=time interval string
%       PRO: Internal function/constant calls are short.
%       CON: create()/parse() function names are not clear.
%     time_interval_string -- IMPLEMENTED
%       PRO: create/parse function names can be short.



  %###################
  %###################
  % PRIVATE CONSTANTS
  %###################
  %###################
  properties(Constant, Access=private)
    % Regular expression for date yyyymmdd (8 digits).
    DATE_RE = '[0-9][0-9][0-9][0-9][01][0-9][0-3][0-9]';

    % Regular expression for date+time yyyymmddThhmmss
    % (8+1+6=15 digits & character)
    %
    % NOTE: DATETIME_RE not same as glob.attr. "Datetime", but is a component
    % of it.
    DATETIME_RE = '[0-9][0-9][0-9]{6,6}T[0-9]{6,6}';
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Create time interval string.
    %
    function timeIntervalStr = create(Dt1, Dt2, timeIntervalFormat)
      irf.dt.assert_UTC(Dt1)
      irf.dt.assert_UTC(Dt2)

      if strcmp(timeIntervalFormat, 'DAY')
        irf.dt.assert_UTC_midnight(Dt1)
        irf.dt.assert_UTC_midnight(Dt2)
        expDt2 = datetime(Dt1) + caldays(1);
        assert(isequal(Dt2, expDt2))

        dateVec1 = datevec(Dt1);

        timeIntervalStr = sprintf(...
          '%04i%02i%02i', ...
          dateVec1(1:3));

      elseif strcmp(timeIntervalFormat, 'DAY_TO_DAY')
        irf.dt.assert_UTC_midnight(Dt1)
        irf.dt.assert_UTC_midnight(Dt2)

        Dt2 = Dt2 - caldays(1);   % Decrease by one day, since end day is inclusive.

        dateVec1 = datevec(Dt1);
        dateVec2 = datevec(Dt2);

        timeIntervalStr = sprintf(...
          '%04i%02i%02i-%04i%02i%02i', ...
          dateVec1(1:3), dateVec2(1:3));

      elseif strcmp(timeIntervalFormat, 'SECOND_TO_SECOND')
        timeIntervalStr = sprintf(...
          '%04i%02i%02iT%02i%02i%02i-%04i%02i%02iT%02i%02i%02i', ...
          datevec(Dt1), datevec(Dt2));

      elseif strcmp(timeIntervalFormat, 'NO_TIME_INTERVAL')
        NAT = irf.dt.UTC('NaT');
        assert(isequaln(Dt1, NAT))
        assert(isequaln(Dt2, NAT))

        timeIntervalStr = [];

      else

        error('Can not interpret timeIntervalFormat="%s".', timeIntervalFormat)
      end

    end



    % Parse time interval string.
    %
    function [Dt1, Dt2, timeIntervalFormat, bSuccess] = parse(timeIntervalStr)
      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
        {solo.adm.dsfn.time_interval_str.DATE_RE}, 'permit non-match');
      if perfectMatch
        Dt1                = solo.adm.dsfn.time_interval_str.day_str_to_DT(subStrCa{1});
        Dt2                = Dt1 + caldays(1);   % Increment by 1 day.
        timeIntervalFormat = 'DAY';
        bSuccess           = true;
        return
      end

      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
        {solo.adm.dsfn.time_interval_str.DATE_RE, '-', ...
        solo.adm.dsfn.time_interval_str.DATE_RE}, 'permit non-match');
      if perfectMatch
        Dt1                = solo.adm.dsfn.time_interval_str.day_str_to_DT(subStrCa{1});
        Dt2                = solo.adm.dsfn.time_interval_str.day_str_to_DT(subStrCa{3});
        Dt2                = Dt2 + caldays(1);   % Increment by 1 day since end day is inclusive.
        timeIntervalFormat = 'DAY_TO_DAY';
        bSuccess           = true;
        return
      end

      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
        {solo.adm.dsfn.time_interval_str.DATETIME_RE, '-', ...
        solo.adm.dsfn.time_interval_str.DATETIME_RE}, 'permit non-match');
      if perfectMatch
        Dt1                = solo.adm.dsfn.time_interval_str.second_str_to_DT(subStrCa{1});
        Dt2                = solo.adm.dsfn.time_interval_str.second_str_to_DT(subStrCa{3});
        timeIntervalFormat = 'SECOND_TO_SECOND';
        bSuccess           = true;
        return
      end

      Dt1                = irf.dt.UTC('NaT');
      Dt2                = irf.dt.UTC('NaT');
      timeIntervalFormat = [];
      bSuccess           = false;
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Utility function
    function Dt = second_str_to_DT(s)
      dateVec = str2double({s(1:4), s(5:6), s(7:8), s(10:11), s(12:13), s(14:15)});
      Dt      = irf.dt.UTC(dateVec);
    end



    % Utility function
    function Dt = day_str_to_DT(s)
      dateVec      = str2double({s(1:4), s(5:6), s(7:8)});
      dateVec(4:6) = [0, 0, 0];
      Dt           = irf.dt.UTC(dateVec);
    end



  end    % methods(Static, Access=private)



end
