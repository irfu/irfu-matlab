function res = caa_append_data(data,app)
%CAA_APPEND_DATA  concatenate two datasets
%
% res = caa_append_data(data,app)
%       append APP to the end of DATA and fill the gap in between
%
% See also CAA_FILL_GAPS
%
% $Id$

% Copyright 2006 Yuri Khotyaintsev

if size(data,2) ~= size(app,2), error('data has a different dimension'), end

res = data;

ii = find(app(:,1)>data(end,1));
if isempty(ii), return, end
app = app(ii,:);

% Fill the gap between the datasets
res = caa_fill_gaps(res,app(1,1));

res = [res; app];
