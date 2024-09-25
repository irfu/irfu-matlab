%
% Parse time interval string. The function supports
% solo.adm.dsfn.DatasetFilename.
%
% NOTE: Tested directly via solo.adm.dsfn.time_interval_str___UTEST and
% indirectly bia solo.adm.dsfn.DatasetFilename___UTEST.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [Dt1, Dt2, timeIntervalFormat] = parse_time_interval_str(timeIntervalStr)
% PROPOSAL: Return special value(s) for illegal string.
%   PRO: Useful for rigorously distinguishing datasets/non-datasets.
%   PROPOSAL: Return struct (not class!).
%   PROPOSAL: timeIntervalFormat == special value, e.g. [].
%
% PROPOSAL: New abbrev. TIFID = Time Interval Format ID
%
% PROPOSAL: Class with static methods for parse_time_interval_str() and
%           create_time_interval_str().
%   PRO: Natural to group parse+create functions.
%   PRO: Natural to have test for class.
%   TODO-DEC: Name?
%     tis=time interval string
%     time_interval_string
%       PROPOSAL: Shorter function names since class name is informative.
%       time_interval_string.create()
%       time_interval_string.parse()

% yyyymmdd (8 digits).
DATE_RE     = '[0-9][0-9][0-9][0-9][01][0-9][0-3][0-9]';

% NOTE: DATETIME_RE not same as glob.attr. Datetime, but a component of it.
% yyyymmddThhmmss (8+1+6=15 digits/T)
DATETIME_RE = '[0-9][0-9][0-9]{6,6}T[0-9]{6,6}';



[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATE_RE}, 'permit non-match');
if perfectMatch
  Dt1 = day_str_to_DT(subStrCa{1});
  Dt2 = Dt1 + caldays(1);   % Increment by 1 day.
  timeIntervalFormat = 'DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATE_RE, '-', DATE_RE}, 'permit non-match');
if perfectMatch
  Dt1 = day_str_to_DT(subStrCa{1});
  Dt2 = day_str_to_DT(subStrCa{3});
  Dt2 = Dt2 + caldays(1);   % Increment by 1 day since end day is inclusive.
  timeIntervalFormat = 'DAY_TO_DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATETIME_RE, '-', DATETIME_RE}, 'permit non-match');
if perfectMatch
  Dt1 = second_str_to_DT(subStrCa{1});
  Dt2 = second_str_to_DT(subStrCa{3});
  timeIntervalFormat = 'SECOND_TO_SECOND';
  return
end

error('Can not interpret time interval string "%s".', timeIntervalStr)
end



% Utility function
function Dt = second_str_to_DT(s)
dateVec = str2double({s(1:4), s(5:6), s(7:8), s(10:11), s(12:13), s(14:15)});
Dt      = datetime(dateVec, 'TimeZone', 'UTCLeapSeconds');
end



% Utility function
function Dt = day_str_to_DT(s)
dateVec      = str2double({s(1:4), s(5:6), s(7:8)});
dateVec(4:6) = [0, 0, 0];
Dt           = datetime(dateVec, 'TimeZone', 'UTCLeapSeconds');
end