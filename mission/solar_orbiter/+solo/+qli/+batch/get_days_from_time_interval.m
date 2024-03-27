%
% Wrapper around solo.qli.offgen.generate_quicklooks() for generating quicklooks
% for all dates within a specified time interval.
%
%
% NOTE: Several arguments are designed to partly handle arguments deriving from
%       the bash/the OS and are therefore on a string format.
%
%
% ARGUMENTS
% =========
% beginDayUtcInclStr, endDayUtcExclStr
%       Strings on format YYYY-MM-DD.
%       Beginning and end of time interval for which quicklooks should be
%       generated.
%       beginDayUtcInclStr is INCLUSIVE. endDayUtcExclStr is EXCLUSIVE, i.e.
%       the day AFTER the last day of the time interval.
%
%
% RETURN VALUES
% =============
% (None)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function DaysDtArray = get_days_from_time_interval(...
  beginDayUtcInclStr, endDayUtcExclStr)

assert(ischar(beginDayUtcInclStr))
assert(ischar(endDayUtcExclStr))

BeginDayInclDt = solo.qli.utils.umdt(beginDayUtcInclStr);
EndDayExclDt   = solo.qli.utils.umdt(endDayUtcExclStr);

% NOTE: Indirectly assertion on the string timestamps.
solo.qli.utils.assert_UTC_midnight_datetime(BeginDayInclDt)
solo.qli.utils.assert_UTC_midnight_datetime(EndDayExclDt)



% IMPLEMENTATION NOTE: Subtracting one day from argument for end timestamp to
% ensure that it is in accordance with the definition of the corresponding
% argument.

% IMPLEMENTATION NOTE: Needs to use "caldays()" not "days()" for handling leap
%                      seconds.
EndDayInclDt = EndDayExclDt - caldays(1);



% Construct array of timestamps, where every timestamp represents one day
% beginning on that timestamp.
% -----------------------------------------------------------------------
% IMPLEMENTATION NOTE: Needs to use "caldays()" not "days()" for handling leap
%                      seconds.
DaysDtArray = [BeginDayInclDt:caldays(1):EndDayInclDt]';

end
