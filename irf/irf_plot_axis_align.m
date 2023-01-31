function irf_plot_axis_align(nsubplot,h)
%IRF_PLOT_AXIS_ALIGN  align axis of different subplots
% usefull when due to colorbars and other stuff axis gets disaligned
%
% irf_plot_axis_align                   align in current figure all handles
% irf_plot_axis_align(handles)          align only handles
% irf_plot_axis_align(nsbuplot)         align with respect to subplot nsubplot
% irf_plot_axis_align(nsbuplot,handles)
%
% Input:
%      handles: axes handles which to align, if not given, align all axes
%     nsbuplot: subplot with respect to which align axis
%               if none given align to fit all
%

narginchk(0,2)

% calculate the size of axis to which align
if nargin==0 % no subplot number given, find smallest common limit for all axis
  h=irf_plot_get_subplot_handles;
  hpostmp=zeros(numel(h),4);
  for j=1:length(h)
    hpostmp(j,:)=get(h(j),'position');
  end
  xpos=[max(hpostmp(:,1)) min(hpostmp(:,3))];
elseif nargin==1 % subplot number or handles to align are specified
  if any(ishandle(nsubplot)) % syntax irf_plot_axis_align(handles) has been used
    h=nsubplot;
    hpostmp=zeros(numel(h),4);
    for j=1:length(h)
      hpostmp(j,:)=get(h(j),'position');
    end
    xpos=[max(hpostmp(:,1)) min(hpostmp(:,3))];
  else % syntax irf_plot_axis_align(nsubplot) has been used
    h=irf_plot_get_subplot_handles;
    hpostmp=get(h(nsubplot),'position');
    xpos=hpostmp([1 3]);
  end
else
  hpostmp=get(h(nsubplot),'position');
  xpos=hpostmp([1 3]);
end

% align subplots
for jj=1:length(h)
  hpos=get(h(jj),'position');
  hpos(1)=xpos(1);
  hpos(3)=xpos(2);
  set(h(jj),'position',hpos);
end

%align ylabels
irf_plot_ylabels_align(h)
