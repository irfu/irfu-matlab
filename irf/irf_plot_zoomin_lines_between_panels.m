function [hline1,hline2]=irf_plot_zoomin_lines_between_panels(h1,h2)
% IRF_PLOT_ZOOMIN_LINES_BETWEEN_PANELS connect with zoomin lines two panels
%
% IRF_PLOT_ZOOMIN_LINES_BETWEEN_PANELS(H1,H2)
%	Connects the upper corners of panel with axes handle H2 with lines
%	to the corresponding time on the bottom of panel with axes handle H1.
%	Thus clearly showing to which time interval panel H2 corresponds in panel H1.
%
% [hline1,hline2]=IRF_PLOT_ZOOMIN_LINES_BETWEEN_PANELS(H1,H2) 
%	Returns handles to both lines

if nargin == 0 % show only help
    help irf_plot_zoomin_lines_between_panels;
    return
elseif nargin ~= 2 % number of input argument is not 2
	irf_log('fcal','WARNING: INCORRECT NUMBER OF ARGUMENTS!');
	return;
end

hh2_xlim=get(h2,'xlim');
hh2_ylim=get(h2,'ylim');
hh1_ylim=get(h1,'ylim');
[zoom_nfu_hh1_x1,zoom_nfu_hh1_y]= ds2nfu(h1,hh2_xlim(1),hh1_ylim(1));
[zoom_nfu_hh1_x2,zoom_nfu_hh1_y]= ds2nfu(h1,hh2_xlim(2),hh1_ylim(1));
hh2_pos=get(h2,'position');
hline1=annotation('line',[hh2_pos(1) zoom_nfu_hh1_x1],[hh2_pos(2)+hh2_pos(4) zoom_nfu_hh1_y]);
hline2=annotation('line',[hh2_pos(1)+hh2_pos(3) zoom_nfu_hh1_x2],[hh2_pos(2)+hh2_pos(4) zoom_nfu_hh1_y]);

if nargout==0
	clear hline1 hline2;
end