function [new_data,offset] = caa_corof_adc(data,start_time,dt,whip)
%caa_corof_adc correct the ADC offset in data by removing the average
%   [new_data,offset] = caa_corof_adc(data) 
%   [new_data,offset] = caa_corof_adc(data,whip) 
%   [new_data,offset] = caa_corof_adc(data,start_time,dt) 
%   [new_data,offset] = caa_corof_adc(data,start_time,dt,whip) 
%
%   $Revision$  $Date$
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if data(end,1)-data(1,1)<60 %interval is shorter then 1 min
	irf_log('proc','Time interval too short, processing may give strange results...')
end

new_data = data;
if nargin==3 | nargin==4
	test_data = irf_tlim(data,start_time + [0 dt]);
	if isempty(test_data)
		irf_log('proc','ADC offset: defaulting to the whole data interval')
		test_data = data;
	end
else
	test_data = data;
end

if nargin==2, whip = start_time; end
if nargin==2 | nargin==4
	irf_log('proc','ADC offset: not using times with Whisper pulses')
	test_data = caa_rm_blankt(test_data, whip); 
end

test_data = test_data(:,2);

ii = find(~isnan(test_data));
if isempty(ii), offset = [];
else, 
	offset = mean(test_data(ii)); clear test_data
	new_data(:,2) = data(:,2) - offset;
end

