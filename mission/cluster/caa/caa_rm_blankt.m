function res = caa_rm_blankt(data,time_int,mode)
%CAA_RM_BLANKT  remove time intervals from data (sets to NaN)
%
% res = caa_rm_blankt(data,time_int,[mode])
%
% Input:
%   data     - Cluster AV format
%   time_int - time intervals (ISDAT epoch)
%   mode     - if 1 remove data OUTSIDE the time intervals
%
% Example:
%   P_NO_WHI = caa_rm_blankt(P1,WHIP1);
%   E_WAKE = caa_rm_blankt(diEs1p34,LOWAKE1p34);
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,3)

res = data;
if isempty(data) || isempty(time_int), return, end

time_int = sort(time_int,1);

for j=1:size(time_int,1)
  if nargin == 2 || mode ~= 1
    res( data(:,1)>=time_int(j,1) & data(:,1)<=time_int(j,2) ,2:end) = NaN;
  else
    if j==1
      res( data(:,1)<time_int(j,1) ,2:end) = NaN;
    else
      res( data(:,1)<time_int(j,1) & data(:,1)>time_int(j-1,2) ,...
        2:end) = NaN;
    end
    if j==size(time_int,1)
      res( data(:,1)>time_int(j,2) ,2:end) = NaN;
    end
  end
end
