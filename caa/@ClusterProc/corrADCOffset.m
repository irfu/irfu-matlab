function [new_data,offset] = corrADCOffset(cp,data)
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
offset = mean(data(:,2));
new_data(:,2) = data(:,2) - offset;
