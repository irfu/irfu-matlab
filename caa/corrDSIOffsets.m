function new_data = corrDSIOffsets(data,Dx,Dy,Da)
%corrDSIOffsets correct offsets in DSI
% newdata = Da*(data_x-Dx, data_y-Dy)
%
% $Id$
%

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)



new_data = data;
new_data(:,2) = data(:,2) - Dx;
new_data(:,3) = data(:,3) - Dy;
new_data(:,2:3) = Da*new_data(:,2:3);
