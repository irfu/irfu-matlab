%
% Log standardized execution speed test. Logs 
% (1) wall time between calls to tic() (called before this function), and toc()
%     (called by this function).
% (2) (optionally) wall time per arbitrary unit of "something".
%
%
% ARGUMENTS
% =========
% L
%       bicas.Logger object.
% codeName
%       Name of code being speed tested, e.g. function name, informal name of a
%       segment of code.
% tTicToc
%       Value returned from "tic()".
% --- Optional: Both arguments or neither. ---
% nUnits
%       Size of the amount of processing done in the form of repeated
%       "somethings" done.
%       Ex: number of CDF records.
% unitName
%       Singular (grammatically). E.g. "record" (not "records").
% --
% NOTE: Function does not follow convention that L is the last argument.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-27.
%
function log_speed_profiling(L, codeName, tTicToc, nUnits, unitName)
% PROPOSAL: Permit multiple units.
%   Ex: Downsampled: Nbr of DSR CDF records (bins), nbr of OSR CDF records.
%   PROPOSAL: Cell array
%       {iUnit}{1} = nUnits
%       {iUnit}{2} = unitName
%       TODO-DEC: Log format? Much data together with long row prefixes.
%         Ex: 2021-05-20
%         2021-05-20T18:10:08 -- DEBUG -- SPEED -- bicas.proc.dsr.get_downsampling_bins: 0.253965 [s], 7.23335e-07 [s/OSR record], 351103 [OSR records] (wall time)
%         2021-05-20T18:10:08 -- DEBUG -- SPEED -- bicas.proc.dsr.get_downsampling_bins: 0.254625 [s], 0.000705332 [s/DSR record], 361 [DSR records] (wall time)
%         2021-05-20T18:10:09 -- DEBUG -- SPEED -- bicas.proc.L2L3.process_L2_to_L3: 0.461114 [s], 1.31333e-06 [s/OSR record], 351103 [OSR records] (wall time)
%         2021-05-20T18:10:09 -- DEBUG -- SPEED -- bicas.proc.L2L3.process_L2_to_L3: 0.461797 [s], 0.00127922 [s/DSR record], 361 [DSR records] (wall time)
%   PROPOSAL: One row per unit, repeat everything else.
%         bicas.proc.dsr.get_downsampling_bins: 0.253965 [s], 7.23335e-07 [s/OSR record], 351103 [OSR records] (wall time)
%         bicas.proc.dsr.get_downsampling_bins: 0.254625 [s], 0.000705332 [s/DSR record], 361 [DSR records] (wall time)
%       PRO: Rows searchable (grep).
%       NOTE: Crude, but still improvement over status quo.

    wallTimeSec = toc(tTicToc);

    if nargin == 3
        infoStr = sprintf('%s: %g [s]', codeName, wallTimeSec);
    elseif nargin == 5
        % NOTE: Adds "s" after unit to get plural.
        % IMPLEMENTATION NOTE: nUnits might be an integer. Must convert to
        % double for division to work.
        %   Ex: bicas.proc.dsr.get_downsampling_bins(): nBins
        infoStr = sprintf('%s: %g [s], %g [s/%s], %g [%ss]', ...
            codeName, ...
            wallTimeSec, wallTimeSec/double(nUnits), unitName, nUnits, unitName);
    else
        error('BICAS:Assertion', 'Illegal number of arguments.')
    end
    
    L.logf('debug', ...
        'SPEED -- %s (wall time)', ...
        infoStr)
end
