%
% Generate quicklooks, but specifying method for how to select dates.
%
%
% NOTE: See README.TXT for information on this package.
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.generate_quicklooks_all_types() instead.
%
% NOTE: This function is designed to be called from bash/the OS.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.generate_quicklooks_all_types() instead.
%
%
% IMPLEMENTATION NOTES
% ====================
% This function effectively wraps two different MATLAB functions so that one can
% write just ONE non-trivial bash script which calls/wraps this function instead
% of multiple. Such a bash script may add "non-trivial" functionality which one
% may not want to implement in multiple bash wrapper scripts, e.g. (1) copying
% output files to an intermediate directory and/or (2) logs the output.
% --
% The functionality (switch-case, shuffling arguments) could in principle have
% been implemented in such a bash wrapper script, but it is better to minimize
% the code in bash wrapper scripts.
%
%
% ARGUMENTS
% =========
% modeId
%     String constant which determines which other MATLAB function to delegate
%     to.
% varargin
%     Arguments passed on to other MATLAB function. See implementation.
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_bash_any(modeId, varargin)
% PROPOSAL: Option for returning help text.
%   PRO: Eliminates duplication of documentation in bash wrapper script.
%     CON: Help text in this function would duplicate documentation in the two
%          functions called.
%   CON: Bash wrapper script needs to be aware of syntax for generating help
%        text so that it does not log or create intermediate directory etc.

switch(modeId)
  case 'TIME_INTERVAL'
    solo.qli.offgen.generate_quicklooks_bash_time_interval(varargin{:})

  case 'GENERATE_FROM_LOGS'
    solo.qli.offgen.generate_quicklooks_bash_from_logs(varargin{:})

  otherwise
    error('Illegal modeId="%s"', modeId)
end

end
