%
% Log standardized execution speed test. Logs wall time between calls to tic() and toc().
%
%
% ARGUMENTS
% =========
% L        : bicas.logger object.
% codeName : Name of code being speed tested, e.g. function name.
% tTicToc  : Value returned from "tic".
% --
% (Optional: Both arguments or neither.)
% nUnits   : E.g. number of records.
% unitName : Singular. E.g. "record".
%
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-27.
%
function log_speed_profiling(L, codeName, tTicToc, nUnits, unitName)
    %
    % PROPOSAL: Log method for logging speed tests.
    %   TODO-DEC: Log message format?
    %       -- INFO -- Speed test (wall time): <code name>: x [s]
    %       -- INFO -- Speed test (wall time): <code name>: x [s], y [s/<unit>]
    %   TODO-DEC: Name
    %       ~log
    %       ~speed
    %       ~test
    %       ~speed result
    %       ~profiling
    %
    % NOTE: Function does not follow convention that L is the last argument.
    
    wallTimeSec = toc(tTicToc);
    
    if nargin == 3
        infoStr = sprintf('%s: %g [s]', codeName, wallTimeSec);
    elseif nargin == 5
        % NOTE: Adds "s" after unit to get plural.
        infoStr = sprintf('%s: %g [s], %g [s/%s], %g [%ss]', codeName, wallTimeSec, wallTimeSec/nUnits, unitName, nUnits, unitName);
    else
        error('Illegal number of arguments.')
    end
    
    L.logf('debug', ...
        'SPEED -- %s (wall time)', ...
        infoStr)
end
