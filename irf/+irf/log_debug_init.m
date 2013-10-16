function f = log_debug_init(logLevel)
%LOG_DEBUG_INIT  initialize debug logging
%
%
% debug_log = irf.log_debug_init();
%   returns function debug_log(msg) with default logging settings
%
% debug_log = irf.log_debug_init(level);
%   returns function debug_log(msg) with local logging settings
%   level = 0 - disabled, level = 1 - enabled
%
% Example:
%   debug_log = irf.log_debug_init(1);
%   debug_log('My debug message')
%
% See also : IRF.LOG

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin == 1
  localLogLevel = logLevel;
elseif nargin == 0
  localLogLevel = irf.log();
  if localLogLevel<4, localLogLevel = 0; end
else
    error('incorrect number of input parameters')
end

if localLogLevel
  f = @(x) irf.log(4,x);
else
  f = @(x) [];
end
  