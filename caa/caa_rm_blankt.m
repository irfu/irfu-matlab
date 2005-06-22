function res = caa_rm_blankt(data,time_int)
%caa_rm_blankt remove time intervals from data (sets to NaN)
% res = caa_rm_blankt(data,time_int)
% Input:
% data - Cluster AV format
% time_int - time intervals (ISDAT epoch)
% 
% Example:
% P_NO_WHI=caa_rm_blankt(P1,WHIP1);
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(2,2,nargin))

res = data;
if isempty(data), return, end

for j=1:size(time_int,1)
	res(find(data(:,1)>=time_int(j,1) & data(:,1)<=time_int(j,2)),2:end) = NaN;
end
