%
% Log standardized execution speed test. Logs wall time between calls to tic()
% (called before this function, and toc() (called by this function).
%
%
% ARGUMENTS
% =========
% L        : bicas.logger object.
% codeName : Name of code being speed tested, e.g. function name.
% tTicToc  : Value returned from "tic()".
% --
% --- Optional: Both arguments or neither. ---
% nUnits   : E.g. number of records.
% unitName : Singular. E.g. "record".
%
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-27.
%
function log_speed_profiling(L, codeName, tTicToc, nUnits, unitName)

    % NOTE: Function does not follow convention that L is the last argument.

    wallTimeSec = toc(tTicToc);
    
    if nargin == 3
        infoStr = sprintf('%s: %g [s]', codeName, wallTimeSec);
    elseif nargin == 5
        % NOTE: Adds "s" after unit to get plural.
        % IMPLEMENTATION NOTE: nUnits might be an integer. Must convert to
        % double for division to work.
        %   Ex: bicas.proc.dwns.get_downsampling_bins(): nBins
        infoStr = sprintf('%s: %g [s], %g [s/%s], %g [%ss]', ...
            codeName, ...
            wallTimeSec, wallTimeSec/double(nUnits), unitName, nUnits, unitName);
    else
        error('Illegal number of arguments.')
    end
    
    L.logf('debug', ...
        'SPEED -- %s (wall time)', ...
        infoStr)
end
