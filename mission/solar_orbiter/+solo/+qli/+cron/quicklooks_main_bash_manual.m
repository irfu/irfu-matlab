%
% Wrapper around solo.qli.cron.quicklooks_cron() intended for being used by
% being called from system scripts (e.g. bash) for the purpose of cron jobs and
% manual official processing on brain/spis. The arguments have also been
% designed for this purpose and are therefore all strings.
%
% NOTE: This function is intended to be called from bash/the OS, hence "bash" in
%       the name.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.generate_quicklooks_all_types() instead.
%
%
% ARGUMENTS
% =========
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       One-character strings. Whether to generate ("1" or "0") non-weekly (2h,
%       6h, 24h) quicklooks and/or weekly quicklooks.
%       NOTE: STRINGS.
% beginDayUtcInclStr, endDayUtcExclStr
%       Strings on format YYYY-MM-DD.
%       Beginning and end of time interval for which quicklooks should be
%       generated.
%       beginDayUtcInclStr is INCLUSIVE. endDayUtcExclStr is EXCLUSIVE, i.e.
%       the day AFTER the last day of the time interval.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function quicklooks_main_bash_manual(...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  beginDayUtcInclStr, endDayUtcExclStr)

generateNonweeklyQuicklooks = interpret_argument_flag(generateNonweeklyQuicklooks);
generateWeeklyQuicklooks    = interpret_argument_flag(generateWeeklyQuicklooks);
assert(ischar(beginDayUtcInclStr))
assert(ischar(endDayUtcExclStr))

BeginDayInclDt          = datetime(beginDayUtcInclStr);
EndDayExclDt            = datetime(endDayUtcExclStr);
BeginDayInclDt.TimeZone = 'UTCLeapSeconds';
EndDayExclDt.TimeZone   = 'UTCLeapSeconds';

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



%=====================
% Generate quicklooks
%=====================
solo.qli.cron.quicklooks_main(...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
end



% Interpret argument for main function interface. Intended accept and normalize
% arguments which are either
% (1) MATLAB-friendly (numeric/logical), or
% (2) bash script-friendly (strings).
%
function value = interpret_argument_flag(arg)
assert(isscalar(arg), 'Flag argument is not scalar.')

if     ischar(arg) && arg=='0'
  value = false;
elseif ischar(arg) && arg=='1'
  value = true;
else
  error('Can not interpret argument flag. Illegal format.')
end
end
