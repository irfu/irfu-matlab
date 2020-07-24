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
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-05-31.
%
function stdout_print(msgStr)
    % PROPOSAL: Change name to something analogous with logf.
    
    printStr = EJ_library.str.add_prefix_on_every_row(msgStr, bicas.constants.STDOUT_PREFIX_TBW);
    
    fwrite(1, printStr);    % NOTE: Must print using function that reacts to trailing line feed.
end
