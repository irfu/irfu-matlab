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
% operationId
%       String constant which specifies what to do with list of dates.
% dateSelectionAlgorithmId
%       String constant which specifies which algorithm/method should be used
%       for obtaining dates.
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
  operationId, dateSelectionAlgorithmId, varargin)

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
%
% PROPOSAL: More consistent string constants.
%   TODO-DEC: How/whether to pluralize "LOG" and "FMD"?
%   PROPOSAL:
%     GENERATE_FROM_*: *_TIME_INTERVAL, *_LOGS, *_FMDS
% PROPOSAL: Additional string constant for whether to list dates or generate
%           quicklooks. -- IMPLEMENTED
%   <operationModId>           : LIST, GENERATE
%   <dateSelectionAlgorithmId> : TIME_INTERVAL, LOGS, FMDS


generateNonweeklyQuicklooks = solo.qli.batch.interface.interpret_boolean_flag(generateNonweeklyQuicklooks);
generateWeeklyQuicklooks    = solo.qli.batch.interface.interpret_boolean_flag(generateWeeklyQuicklooks);

DaysDtArray = solo.qli.batch.interface.get_days_from_selected_algorithm(...
  Settings, outputDir, dateSelectionAlgorithmId, varargin);

switch(operationId)
  case 'LIST'
    %============
    % List dates
    %============
    % NOTE: Does not use irf.log().

    % Compensate for that previous log message appears to not end with line
    % feed.
    fprintf('\n')

    nDates = numel(DaysDtArray);
    for i = 1:nDates
      Dt = DaysDtArray(i);
      fprintf('%04i-%02i-%02i\n', Dt.Year, Dt.Month, Dt.Day)
    end
    fprintf('Number of dates: %i\n', nDates)

  case 'GENERATE'
    %=====================
    % Generate quicklooks
    %=====================
    solo.qli.batch.generate_quicklooks(...
      Settings.irfLogoPath, Settings.vhtDir, outputDir, ...
      generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
      DaysDtArray, Settings.Gql)

  otherwise
    error('Illegal operationId="%s"', operationId)
end

end
