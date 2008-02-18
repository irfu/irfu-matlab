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
% Dx,Dy can also be time dependant (AV_CLUSTER format)
%
% See also C_EFW_DSI_OFF
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(3,4,nargin))

if length(Dx) == 1
	if nargin ==3
		Da = Dy;
		Dy = imag(Dx);
		Dx = real(Dx);
	end
	new_data = corr_local(data,Dx,Dy,Da);
	return
end

if nargin ==3
	Da = Dy;
	Dy = Dx; Dy(:,2) = imag(Dx(:,2));
	Dx(:,2) = real(Dx(:,2));
elseif any(size(Dx)~=size(Dy))
	error('Dx and Dy must be of the same size')
end

new_data = data;
for in=1:size(Dx,1)
	if in ==1, ii = find( data(:,1)<Dx(2,1) );
	elseif in == size(Dx,1), ii = find( data(:,1)>=Dx(end,1) );
	else ii = find( data(:,1)<Dx(in+1,1) & data(:,1)>=Dx(in,1) );
	end
	new_data(ii,:) = corr_local(data(ii,:),Dx(in,2),Dy(in,2),Da);
end

function new_data = corr_local(data,Dx,Dy,Da)
new_data = data;
new_data(:,2) = data(:,2) - Dx;
new_data(:,3) = data(:,3) - Dy;
new_data(:,2:3) = Da*new_data(:,2:3);
