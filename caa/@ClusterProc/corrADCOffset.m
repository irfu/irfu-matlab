function [new_data,offset] = corrADCOffset(cp,data,start_time,dt)
%corrADCOffset correct the ADC offset
%   [new_data,offset] = corrADCOffset(data) corrects the ADC offset 
%	in p12 and p34 fignals by removing the average
%
%   $Revision$  $Date$
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if data(end,1)-data(1,1)<60 %interval is shorter then 1 min
	warning('Time interval too short, processing may give strange results...')
end

new_data = data;
if nargin==4
	test_data = av_t_lim(data,start_time + [0 dt]);
	if isempty(test_data)
		test_data = data(:,2);
	else
		test_data = test_data(:,2);
	end
else
	test_data = data(:,2);
end

offset = mean(test_data); clear test_data
new_data(:,2) = data(:,2) - offset;
