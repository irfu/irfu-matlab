function irf_plot_axis_align(nsubplot)
%IRF_PLOT_AXIS_ALIGN  align all axis with respect to given subplot 
%
% usefull when due to colorbars and other stuff axis gets disaligned
%
% irf_plot_axis_align(nsbuplot)
%
% Input:
%     nsbuplot: subplot with respect to which align axis
%
% $Id$

h=irf_plot_get_subplot_handles;
hpostmp=get(h(nsubplot),'position');
xpos=hpostmp([1 3]);
for jj=1:length(h),
    hpos=get(h(jj),'position');
    hpos(1)=xpos(1);
    hpos(3)=xpos(2);
    set(h(jj),'position',hpos);
end
