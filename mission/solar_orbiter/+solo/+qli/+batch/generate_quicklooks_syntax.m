%
% Batch-generate quicklooks, but the caller specifies the algorithm (DASA) for
% how to select dates, e.g. specify range, derive from logs, derive from file
% modification dates etc.
%
% NOTE: This function is primarily intended to be called from
%       solo.qli.batch.generate_quicklooks_shell().
%
%
% ARGUMENTS
% =========
% NOTE: Also see wrapper solo.qli.batch.generate_quicklooks_shell()
%       for description of arguments.
% --
% Config
%       Struct with misc. (mostly) system-dependent configuration.
% Gql
%       See solo.qli.batch.generate_quicklooks().
% Fsr
%       See solo.qli.batch.interface.get_days_from_DASA().
% outputDir
% generateNonweeklyQuicklooks, generateWeeklyQuicklooks
% operationId
%       See solo.qli.batch.generate_quicklooks_shell().
% dasaid
%       See solo.qli.batch.generate_quicklooks_shell().
% dasaArgumentsCa
%       See solo.qli.batch.generate_quicklooks_shell().
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_syntax(...
  Config, Gql, Fsr, outputDir, ...
  generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
  operationId, dasaid, dasaArgumentsCa)

% PROPOSAL: Better name.
%   Is not meant to be called from bash, but almost. The exception is settings.
%   ~bash
%   ~interface
%     CON: Is not true "outer interface" (called from non-MATLAB).
%   ~main
%   ~syntax
%   ~DASA
%   ~DASA arguments
%
% ~PROBLEM/UGLY: Specifying generation-specific arguments also when no
%                QLIs should be generated.

% NOTE: soloDbDirPath is not used by this function.
irf.assert.struct(Config, ...
  {'irfLogoPath', 'vhtDir', 'LogFileDirPatternDict', 'fmdQliDir', 'datasetDirsCa'}, ...
  {'soloDbDirPath'})
assert(islogical(generateNonweeklyQuicklooks) & isscalar(generateNonweeklyQuicklooks))
assert(islogical(generateWeeklyQuicklooks)    & isscalar(generateWeeklyQuicklooks))



irf.log('n', sprintf('operationId = "%s"', operationId))
irf.log('n', sprintf('dasaid      = "%s"', dasaid))



UmdDtArray = solo.qli.batch.interface.get_days_from_DASA(...
  Config.datasetDirsCa, ...
  Config.LogFileDirPatternDict, ...
  Fsr, ...
  Config.fmdQliDir, ...
  dasaid, ...
  dasaArgumentsCa);

switch(operationId)

  case 'LIST'
    %=============================================
    % List dates (do not generate any quicklooks)
    %=============================================
    list_operation(UmdDtArray)

  case 'GENERATE'
    %=====================
    % Generate quicklooks
    %=====================
    solo.qli.batch.generate_quicklooks(...
      Config.irfLogoPath, Config.vhtDir, outputDir, ...
      generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
      UmdDtArray, Gql)

  otherwise
    error('Illegal operationId="%s"', operationId)

end

end



% NOTE: Does not use irf.log().
function list_operation(UmdDtArray)

% Compensate for that preceding log message appears to not end with line feed.
fprintf('\n')

nDates = numel(UmdDtArray);
for i = 1:nDates
  Dt = UmdDtArray(i);
  fprintf('%04i-%02i-%02i\n', Dt.Year, Dt.Month, Dt.Day)
end
fprintf('Number of dates: %i\n', nDates)

end
