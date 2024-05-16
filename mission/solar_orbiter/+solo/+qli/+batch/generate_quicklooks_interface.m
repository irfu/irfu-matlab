%
% Batch-generate quicklooks, but the caller specifies the method for how to
% select dates, e.g. specify range, derive from logs, derive from file
% modification dates etc.
%
% NOTE: This function is intended for batch generating quicklooks, in particular
%       by being called from another, custom MATLAB function which specifies
%       (hardcodes) the relevant system setup and initializes "SolO DB", and
%       which is in turn called by a bash script.
%
%
% ARGUMENTS
% =========
% Settings
%       Struct with misc. (mostly) system-dependent settings.
%       .fmdQliDir
%           Path to already generated quicklooks which can be used for reading
%           at file modification dates. Not the same as output directory.
%           Should be enough to point to the 24h quicklooks.
%       .irfLogoPath
%           See solo.qli.batch.generate_quicklooks().
%       .vhtDir
%           See solo.qli.batch.generate_quicklooks().
%       .LogFileDirPatternDict
%           Dictionary of log file paths with globbing.
%           See solo.qli.batch.interface.get_days_from_logs().
%       .datasetDirsCa
%           Cell array of paths to directories with datasets. Used for FMDS.
%           See solo.qli.batch.interface.get_days_from_DMRQ().
%       .Gql
%           See solo.qli.batch.generate_quicklooks().
%       .Fsr
%           See solo.qli.batch.interface.get_days_from_DASA().
% outputDir
%       Path to output directory.
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
%       One-character strings. Whether to generate ("1" or "0") non-weekly (2h,
%       6h, 24h) quicklooks and/or weekly quicklooks.
%       NOTE: STRINGS.
% operationId
%       String constant which specifies what to do with list of dates.
% dasaid
%       String constant which specifies which DASA (algorithm) should be used
%       for obtaining dates.
% varargin
%       Arguments associated with the selected DASA. See implementation.
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
  operationId, dasaid, varargin)

% PROPOSAL: Better name.
%   Is not meant to be called from bash, but almost. The exception is settings.
%   ~bash
%   ~interface
%     CON: Is not true "outer interface" (called from non-MATLAB).
%   ~main
%   ~syntax
%   ~DASA
% PROPOSAL: Option for returning help text.
%   PRO: Eliminates duplication of documentation in bash wrapper script.
%     CON: Help text in this function would duplicate documentation in the two
%          functions called.
%   CON: Bash wrapper script needs to be aware of syntax for generating help
%        text so that it does not log or create intermediate directory etc.
% PROPOSAL: Omit "objects" (Fsr, Gql) from Settings.
%
% PROPOSAL: Command for generating all quicklooks older than certain FMD.
%
% PROPOSAL: Arguments which are passed on from bash wrapper for indirectly
%           specifying dates etc., and which are not setup-dependent, should be
%           a well-defined set stored in a 1D cell array.
%   CON: Those arguments are not a well-defined set.
%     CON: operationId-dependent arguments passed on to
%          solo.qli.batch.interface.get_days_from_DASA() are a
%          well-defined set.
%   PRO: Makes it more clear which arguments are passed on (copied) to where.
%
% ~PROBLEM/UGLY: Specifying generation-specific arguments also when no
%                quicklooks should be generated.

irf.assert.struct(Settings, ...
  {'irfLogoPath', 'vhtDir', 'LogFileDirPatternDict', 'fmdQliDir', 'datasetDirsCa', 'Gql', 'Fsr'}, {})



generateNonweeklyQuicklooks = solo.qli.batch.interface.interpret_boolean_flag(generateNonweeklyQuicklooks);
generateWeeklyQuicklooks    = solo.qli.batch.interface.interpret_boolean_flag(generateWeeklyQuicklooks);

irf.log('n', sprintf('operationId = "%s"', operationId))
irf.log('n', sprintf('dasaid      = "%s"', dasaid))



DaysDtArray = solo.qli.batch.interface.get_days_from_DASA(...
  Settings.datasetDirsCa, ...
  Settings.LogFileDirPatternDict, ...
  Settings.Fsr, ...
  Settings.fmdQliDir, ...
  dasaid, ...
  varargin(:));

switch(operationId)

  case 'LIST'
    %=============================================
    % List dates (do not generate any quicklooks)
    %=============================================
    list_operation(DaysDtArray)

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



% NOTE: Does not use irf.log().
function list_operation(DaysDtArray)

% Compensate for that preceding log message appears to not end with line feed.
fprintf('\n')

nDates = numel(DaysDtArray);
for i = 1:nDates
  Dt = DaysDtArray(i);
  fprintf('%04i-%02i-%02i\n', Dt.Year, Dt.Month, Dt.Day)
end
fprintf('Number of dates: %i\n', nDates)

end