function data = caa_append_data(data,app)
%CAA_APPEND_DATA  concatenate two datasets
%
% res = caa_append_data(data,app)
%       append APP to the end of DATA and fill the gap in between
%
% See also CAA_FILL_GAPS
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,2,nargin))

if size(data,2) ~= size(app,2), error('data has a different dimension'), end

if app(1) <= data(end,1)
	ii = find(app(:,1)>data(end,1));
	if isempty(ii), return, end
	app = app(ii,:);
end

% Fill the gap between the datasets
data = caa_fill_gaps(data,app(1,1));

data = [data; app];
