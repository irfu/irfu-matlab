function res = irf_lstyle(n)
%IRF_LSTYLE    Generate line style
%
% RES = IRF_LSTYLE(N)
%       Generate line style for curve number N
%
% Example:
%       for i=1:5
%         plot(t,x(:,i),irf_lstyle(i)); hold on
%       end
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

colors = {'k' 'r' 'g' 'b' 'm' 'c' 'y'};
markers = {'-' '--' '.-' '-.' 'o-'};

mr = floor(n/length(colors));
cl = n - mr*length(colors);

res = sprintf('%s%s',colors{cl+1},markers{mr+1});
