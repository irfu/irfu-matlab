%
% Generate human-readable "context strings" to put at bottom of plots.
%
% IMPLEMENTATION NOTES
% ====================
% Could be split up into two functions, one per string. Is not since both
% strings can be seen as complementary, or as top and bottom string. Makes it
% more natural if future updates information between the strings.
%
%
% ARGUMENTS
% =========
% soloPosTSeries, earthPosTSeries
%       "Complete" TSeries with data.
% Tint
%       Selected time interval.
%
%
% RETURN VALUES
% =============
% Human-readable string.
%       One-row.
%       Empty string if no data for time interval.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [soloStr, earthStr] = context_info_strings(soloPosTSeries, earthPosTSeries, Tint)
    % TODO-DEC: Function name?
    %       Context, info, string(s)
    % PROPOSAL: No Tint argument. Caller submits already truncated TSeries.
    %   PRO: One fewer arguments.
    %   CON: Caller has to truncate twice.
    %   CON: Caller might truncate differently for different TSeries.

    assert(isa(soloPosTSeries,  'TSeries'))
    assert(isa(earthPosTSeries, 'TSeries'))
    assert(isa(Tint,            'EpochTT'))

    % NOTE: In principle a lot of execution/time just for obtaining a constant,
    %       but the function is not time critical so should not be a problem.
    Units = irf_units;
    AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

    soloPos = soloPosTSeries.tlim(Tint).data;
    if ~isempty(soloPos)
        % NOTE: Doubled whitespaces.
        soloStr = ['SolO: ', ...
            [' R=', sprintf('%.2f',soloPos(1,1)/AU_KM),'Au, '],...
            [' EcLat=',sprintf('%d',round(soloPos(1,3)*180/pi)),'\circ, '],...
            [' EcLon=',sprintf('%d',round(soloPos(1,2)*180/pi)),'\circ'] ...
        ];
    else
        soloStr = '';
    end

    earthPos = earthPosTSeries.tlim(Tint).data;
    if ~isempty(earthPos)
        earthStr = ['Earth: EcLon=',sprintf('%d',round(earthPos(1,2)*180/pi)),'\circ'];
    else
        earthStr = '';
    end
end
