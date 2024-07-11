%
% Parse time interval string. Function supports
% solo.adm.dsfn.parse_dataset_filename().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [dateVec1, dateVec2, timeIntervalFormat] = parse_time_interval_str(timeIntervalStr)
    % PROPOSAL: Automatic test code.

% yyyymmdd (8 digits).
DATE_RE     = '20[0-9][0-9][01][0-9][0-3][0-9]';

% NOTE: DATETIME_RE not same as glob.attr. Datetime, but component of.
% yyyymmddThhmmss (8+1+6=15 digits/T)
DATETIME_RE = '20[0-9]{6,6}T[0-9]{6,6}';



[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATE_RE}, 'permit non-match');
if perfectMatch
  dateVec1 = date_str_2_dateVec(subStrCa{1});
  dateVec2 = datevec(datetime(dateVec1) + caldays(1));   % Increment by 1 day.
  timeIntervalFormat = 'DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATE_RE, '-', DATE_RE}, 'permit non-match');
if perfectMatch
  dateVec1 = date_str_2_dateVec(subStrCa{1});
  dateVec2 = date_str_2_dateVec(subStrCa{3});
  dateVec2 = datevec(datetime(dateVec2) + caldays(1));   % Increment by 1 day.
  timeIntervalFormat = 'DAY_TO_DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(timeIntervalStr, ...
  {DATETIME_RE, '-', DATETIME_RE}, 'permit non-match');
if perfectMatch
  dateVec1 = date_time_str_2_dateVec(subStrCa{1});
  dateVec2 = date_time_str_2_dateVec(subStrCa{3});
  timeIntervalFormat = 'SECOND_TO_SECOND';
  return
end

error('Can not interpret time interval string "%s".', timeIntervalStr)
end



% Utility function
function dateVec = date_time_str_2_dateVec(s)
dateVec = str2double({s(1:4), s(5:6), s(7:8), s(10:11), s(12:13), s(14:15)});

% NOTE: Is not a check on filename, but on implementation. read_token()
% should guarantee that strings can be parsed as numbers.
%assert(~any(isnan(dateVec)))
end



% Utility function
function dateVec = date_str_2_dateVec(s)
dateVec = str2double({s(1:4), s(5:6), s(7:8)});
dateVec(4:6) = [0, 0, 0];

% NOTE: Is not a check on filename, but on implementation. read_token should
% guarantee that strings can be parsed as numbers.
%assert(~any(isnan(dateVec)))
end