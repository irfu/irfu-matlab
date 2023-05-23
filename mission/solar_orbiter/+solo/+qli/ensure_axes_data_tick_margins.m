%
% ~Utility function that removes duplicated code from plot functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function ensure_axes_data_tick_margins(hAxesArray)
    assert(isa(hAxesArray, 'matlab.graphics.axis.Axes'))

    for i = 1:numel(hAxesArray)
        h = hAxesArray(i);

        h.YLim = solo.qli.ensure_data_tick_margins(...
            [min(h.YTick), max(h.YTick) ], ...
            [    h.YLim(1),    h.YLim(2)] ...
        );
    end
end
