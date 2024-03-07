%
% Generate human-readable "context strings" to put at bottom of quicklooks.
%
%
% IMPLEMENTATION NOTES
% ====================
% Could be split up into two functions, one per string. It has not been split up
% both strings can be seen as complementary, or as a top and bottom string.
% --
% One can use "R=" etc, but then the text comes too close to the UTC date
% (reduce font size?).
%
%
% ARGUMENTS
% =========
% soloPosTSeries, earthPosTSeries
%       TSeries with positions for SolO and Earth.
% Tint
%       Selected time interval.
%
%
% RETURN VALUES
% =============
% Two human-readable one-row (no line feed) strings.
%       Empty string(s) if no corresponding data for time interval.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [soloStr, earthStr] = context_info_strings(soloPosTSeries, earthPosTSeries, Tint)
% PROPOSAL: Move to solo.qli.utils.
%
% PROPOSAL: No Tint argument. Caller submits already truncated TSeries.
%   PRO: One fewer arguments.
%   CON: Caller has to truncate twice.
%   CON: Caller might truncate differently for different TSeries.

assert(isa(soloPosTSeries,  'TSeries'))
assert(isa(earthPosTSeries, 'TSeries'))
assert(isa(Tint,            'EpochTT'))
assert(length(Tint) == 2)

% NOTE: In principle a lot of execution/time just for obtaining a constant,
%       but the function is not time-critical so it should not be a problem.
Units = irf_units;
AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

soloPos = soloPosTSeries.tlim(Tint).data;
if ~isempty(soloPos)
  % NOTE: 2x whitespaces between every value.
  soloStr = sprintf([
    'SolO:', ...
    '  %.2f AU,', ...
    '  EcLat %d\\circ,', ...
    '  EcLon %d\\circ', ...
    ], ...
    soloPos(1,1)/AU_KM, ...
    round(soloPos(1,3)*180/pi), ...
    round(soloPos(1,2)*180/pi) ...
    );
else
  soloStr = '';
end

earthPos = earthPosTSeries.tlim(Tint).data;
if ~isempty(earthPos)
  earthStr = sprintf('Earth:  EcLon %d\\circ', round(earthPos(1,2)*180/pi));
else
  earthStr = '';
end

end
