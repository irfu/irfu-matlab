%
% Print string to MATLAB's stdout.
% 
% NOTE: Can handle strings with many line feeds.
% NOTE: Adds its own constant prefix to row so that wrapper bash script can filter out the rows.
%
%
% ARGUMENTS
% =========
% msgStr : Potentially multi-row string to be printed. NOTE: Multi-row strings must end with line feed.
%       
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31.
%
function stdout_print(msgStr)
% PROPOSAL: Move STDOUT_PREFIX, LOG_PREFIX to ~constants.
%   PRO: Not overridable.
%   PRO: Faster?!
%
% PROPOSAL: Change name to something analogous with logf.

global SETTINGS

STDOUT_PREFIX = SETTINGS.get_fv('STDOUT_PREFIX');



printStr = EJ_library.utils.add_prefix_on_every_row(msgStr, STDOUT_PREFIX);

fwrite(1, printStr);    % NOTE: Must print using function that reacts to trailing line feed.

end
