function [new_data,offset] = corrADCOffset(data,start_time,dt,whip)
%corrADCOffset correct the ADC offset in data by removing the average
%   [new_data,offset] = corrADCOffset(data) 
%   [new_data,offset] = corrADCOffset(data,whip) 
%   [new_data,offset] = corrADCOffset(data,start_time,dt) 
%   [new_data,offset] = corrADCOffset(data,start_time,dt,whip) 
%
%   $Revision$  $Date$
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'caa_corof_adc')

if data(end,1)-data(1,1)<60 %interval is shorter then 1 min
	c_log('proc','Time interval too short, processing may give strange results...')
end

new_data = data;
if nargin==3 | nargin==4
	test_data = av_t_lim(data,start_time + [0 dt]);
	if isempty(test_data)
		c_log('proc','ADC offset: defaulting to the whole data interval')
		test_data = data;
	end
else
	test_data = data;
end

if nargin==2, whip = start_time; end
if nargin==2 | nargin==4
	c_log('proc','ADC offset: not using times with Whisper pulses')
	test_data = blankTimes(test_data, whip); 
end

test_data = test_data(:,2);

offset = mean(test_data(find(~isnan(test_data)))); clear test_data
new_data(:,2) = data(:,2) - offset;
