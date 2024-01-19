%
% Print string to MATLAB's stdout. Adds its own constant prefix to row so that
% wrapper bash script can filter out the rows.
%
% Cf bicas.stdout_printf().
%
% NOTE: Can handle strings with many line feeds.
%
%
% ARGUMENTS
% =========
% msgStr : Potentially multi-row string to be printed.
%          NOTE: Multi-row strings must end with line feed.
%
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-05-31.
%
function stdout_print(msgStr)

printStr = irf.str.add_prefix_on_every_row(...
  msgStr, bicas.const.STDOUT_PREFIX_TBW);

% NOTE: Must print using function that reacts to trailing line feed.
fwrite(1, printStr);
end
