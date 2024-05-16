%
% Batch-generate quicklooks, but the caller specifies the algorithm (DASA) for
% how to select dates, e.g. specify range, derive from logs, derive from file
% modification dates etc.
%
% NOTE: This function is intended for being called from another, custom MATLAB
%       function which specifies (hardcodes) the relevant system setup and
%       initializes "SolO DB", and which is in turn called by a bash script.
%
%
% ARGUMENTS
% =========
% NOTE: Also see wrapper solo.qli.batch.generate_quicklooks_interface_offgen()
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
% dasaid
%       String constant which specifies which DASA (algorithm) should be used
%       for obtaining dates.
% dasaArgumentsCa
%       1D column cell array. Arguments associated with the selected DASA
%       (dasaid).
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



DaysDtArray = solo.qli.batch.interface.get_days_from_DASA(...
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
    list_operation(DaysDtArray)

  case 'GENERATE'
    %=====================
    % Generate quicklooks
    %=====================
    solo.qli.batch.generate_quicklooks(...
      Config.irfLogoPath, Config.vhtDir, outputDir, ...
      generateNonweeklyQuicklooks, generateWeeklyQuicklooks, ...
      DaysDtArray, Gql)

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
