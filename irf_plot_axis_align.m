function irf_plot_axis_align(nsubplot)
%IRF_PLOT_AXIS_ALIGN  align all axis with respect to given subplot
% usefull when due to colorbars and other stuff axis gets disaligned
%
% irf_plot_axis_align(nsbuplot)
%
% Input:
%     nsbuplot: subplot with respect to which align axis
%               if none given align to the panel
%
% $Id$

error(nargchk(0,1,nargin))
h=irf_plot_get_subplot_handles;

% calculate the size of axis to which align
if nargin==0, % no subplot number given
  for j=1:length(h);
    hpostmp(j,:)=get(h(j),'position');
  end
  xpos=[max(hpostmp(:,1)) min(hpostmp(:,3))];
else % subplot number is given
  hpostmp=get(h(nsubplot),'position');
  xpos=hpostmp([1 3]);
end

% align subplots 
for jj=1:length(h),
  hpos=get(h(jj),'position');
  hpos(1)=xpos(1);
  hpos(3)=xpos(2);
  set(h(jj),'position',hpos);
end
