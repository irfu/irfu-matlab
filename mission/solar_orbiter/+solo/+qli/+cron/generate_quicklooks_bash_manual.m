%
% Wrapper around solo.qli.cron.generate_quicklooks() for generating quicklooks
% for all dates within a specified time interval.
%
%
% NOTE: See README.TXT for information on this package.
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.generate_quicklooks_all_types() instead.
%
% NOTE: This function is designed to potentially be called from bash/the OS.
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory.
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
% RETURN VALUES
% =============
% (None)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function generate_quicklooks_bash_manual(...
  outputDir, generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
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
solo.qli.cron.generate_quicklooks(...
  outputDir, generateNonweeklyQuicklooks, generateWeeklyQuicklooks, DaysDtArray)
end



% Interpret argument for main function interface. Intended accept and normalize
% arguments which are either
% (1) MATLAB-friendly (numeric/logical), or
% (2) bash script-friendly (strings).
%
function value = interpret_argument_flag(arg)
% NOTE: num2str() converts string/number-->string.
assert(isscalar(arg), 'Flag argument "%s" is not scalar.',   num2str(arg))
assert(ischar(arg),   'Flag argument "%s" is not a string.', num2str(arg))

if     ischar(arg) && arg=='0'
  value = false;
elseif ischar(arg) && arg=='1'
  value = true;
else
  error('Can not interpret argument flag="%s". Illegal value.', arg)
end
end
