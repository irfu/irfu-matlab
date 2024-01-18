%
% Save summary plot that has already been created as a MATLAB plot figure.
%
% NOTE: Assumes that
%   * figure already exists
%   * filename has already been chosen
% NOTE: No log message
%
% Useful to have this separately so that it can improved separately (?) from the
% slow generation of summary plots.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-08-13.
%
function save_SP_figure(spPath, hFig)
%saveas(hFig, spPath)
%print(spPath, '-r300')

% IMPLEMENTATION NOTE: Need to set some properties to ensure producing the
% same image on both irony and brain (at least set PaperPosition).
%set(hFig, 'PaperPosition', [0,0, 27.94, 72.2683])
set(hFig, 'PaperPosition', [0,0, 40, 40*sqrt(2)])
%get(hFig)   % DEBUG
print(spPath, '-dpng', '-r200')
end
