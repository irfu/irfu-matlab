function plotExy(varargin)
% plotExy(e1,e2,e3...) plots Ex and Ey on two panels.
% supports multiple data inputs
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev
if nargin<1, error('Need input arguments'),end

figure
p1 = subplot(2,1,1);
p2 = subplot(2,1,2);
nd = length(varargin);

for i=1:nd
	e = varargin{i};
	axes(p1)
	plot(e(:,1),e(:,2),getLineStyle(i-1));
	if i==1, hold on, end
	if i==nd
		hold off
		grid
		add_timeaxis
		ylabel('Ex [mV/m]')
	end
	axes(p2)
	plot(e(:,1),e(:,3),getLineStyle(i-1));
	if i==1, hold on, end
	if i==nd
		hold off
		grid
		add_timeaxis
		ylabel('Ey [mV/m]')
	end
end
addPlotInfo
