function res = caa_is_sh_interval(dir_s)
%CAA_IS_SH_INTERVAL  check if directory contains an SH interval
%
% CAA_IS_SH_INTERVAL([DIR_S])
%     Check for .caa_sh_interval in DIR_S on in the current dir (no DIR_S)
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin < 1, dir_s = pwd; end

if exist([dir_s '/.caa_sh_interval'],'file'), res = 1;
else, res = 0;
end