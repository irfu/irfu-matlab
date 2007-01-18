function new_data = caa_corof_dsi(data,Dx,Dy,Da)
%CAA_COROF_DSI  correct offsets in DSI(ISR2)
%
% new_data = caa_corof_dsi(data,Dx,Dy,Da)
% new_data = caa_corof_dsi(data,Ddsi,Damp)
%   Ddsi is complex: Dx = real(Ddsi), Dy = imag(Ddsi)
%
% The followig computation is performed:
% newdata = Da*(data_x-Dx, data_y-Dy)
%
% $Id$

% Copyright 2004-2007 Yuri Khotyaintsev (yuri@irfu.se)

if nargin ==3,
	Da = Dy;
	Dy = imag(Dx);
	Dx = real(Dx);
end

new_data = data;
new_data(:,2) = data(:,2) - Dx;
new_data(:,3) = data(:,3) - Dy;
new_data(:,2:3) = Da*new_data(:,2:3);
