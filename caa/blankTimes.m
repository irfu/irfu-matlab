function res = blankTimes(data,time_int)
%blankTimes remove time intervals from data (sets to NaN)
% res = blankTimes(data,time_int)
% Input:
% data - Cluster AV format
% time_int - time intervals (ISDAT epoch)
% 
% Example:
% P_NO_WHI=blankTimes(P1,WHIP1);
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(2,2,nargin))

res = data;

for j=1:size(time_int,1)
	res(find(data(:,1)>time_int(j,1) & data(:,1)<time_int(j,2)),:) = NaN;
end
