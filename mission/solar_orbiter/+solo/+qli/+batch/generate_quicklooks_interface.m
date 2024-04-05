%
% Generate quicklooks, but the caller specifies the method for how to
% select dates.
%
% NOTE: The function is intended for batch generating quicklooks, in particular
%       by being called from another, custom MATLAB function which specifies
%       (hardcodes) the relevant system setup and initializes "SolO DB", and
%       which is in turn called by a bash script.
%
%
% ARGUMENTS
% =========
% Settings
%       Struct with misc. system-dependent settings.
% outputDir
%       Path to output directory.
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       One-character strings. Whether to generate ("1" or "0") non-weekly (2h,
%       6h, 24h) quicklooks and/or weekly quicklooks.
%       NOTE: STRINGS.
% dateSelectionAlgorithmId
%       String constant which the date selection algorithm/method.
% varargin
%     Arguments associated with the selected date selection algorithm. See
%     implementation.
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_interface(...
  Settings, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  dateSelectionAlgorithmId, varargin)

% PROPOSAL: Better name.
%   Is not meant to be called from bash, but almost. The exception is settings.
%   ~bash
%   ~interface
%   ~main
%   ~syntax
% PROPOSAL: Option for returning help text.
%   PRO: Eliminates duplication of documentation in bash wrapper script.
%     CON: Help text in this function would duplicate documentation in the two
%          functions called.
%   CON: Bash wrapper script needs to be aware of syntax for generating help
%        text so that it does not log or create intermediate directory etc.
% PROPOSAL: More consistent string constants.
%   TODO-DEC: How/whether to pluralize "LOG" and "FMD"?

generateNonweeklyQuicklooks = interpret_argument_flag(generateNonweeklyQuicklooks);
generateWeeklyQuicklooks    = interpret_argument_flag(generateWeeklyQuicklooks);

switch(dateSelectionAlgorithmId)
  case 'TIME_INTERVAL'
    DaysDtArray = solo.qli.batch.get_days_from_time_interval(varargin{:});

  case 'GENERATE_FROM_LOGS'
    DaysDtArray = solo.qli.batch.get_days_from_logs(Settings, varargin{:});

  case 'GENERATE_FROM_FMDS'
    DaysDtArray = solo.qli.batch.get_days_from_FMDs(Settings, varargin{:});

  otherwise
    error('Illegal dateSelectionAlgorithmId="%s"', dateSelectionAlgorithmId)
end



%=====================
% Generate quicklooks
%=====================
solo.qli.batch.generate_quicklooks(...
  Settings.irfLogoPath, Settings.vhtDir, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  DaysDtArray, Settings.Gql)

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
