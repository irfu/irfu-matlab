%
% ARGUMENTS
% =========
%
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
function [soloStr] = context_info_strings(soloPosKmTSeries, Tint)
    % TODO-DEC: Name?
    %       Context, info, string(s)
    % PROPOSAL: No Tint argument.

    % NOTE: In principle a lot of execution/time just for obtaining a constant,
    %       but the function is not time critical so should not be a problem.
    Units = irf_units;
    AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

    soloPosKm = soloPosKmTSeries.tlim(Tint).data;
    if ~isempty(soloPosKm)
        % NOTE: Doubled whitespaces.
        soloStr = ['SolO: ', ...
            [' R=', sprintf('%.2f',soloPosKm(1,1)/AU_KM),'Au, '],...
            [' EcLat=',sprintf('%d',round(soloPosKm(1,3)*180/pi)),'\circ, '],...
            [' EcLon=',sprintf('%d',round(soloPosKm(1,2)*180/pi)),'\circ'] ...
        ];
    else
        soloStr = char();
    end

end
