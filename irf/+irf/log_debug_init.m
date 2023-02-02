function f = log_debug_init(logLevel)
%LOG_DEBUG_INIT  initialize debug logging
%
% debug_log = irf.log_debug_init(flag);
%   returns function handle 'debug_log'
%
%   If flag==false (default) - debugging disabled
%   If flag==true  - debugging enabled
%
%   debug_log(msg) display debug message 'msg'
%
%   Default logging settings from irf.log apply also to debug_log
%
% Example:
%   debug_log = irf.log_debug_init(true);
%   debug_log('My debug message')
%
% See also : IRF.LOG

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin == 0
  logLevel = false;
elseif nargin > 1
  error('incorrect number of input parameters')
end

if logLevel
  f = @(x) irf.log('debug',x);
else
  f = @(x) [];
end
