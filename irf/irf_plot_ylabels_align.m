function irf_plot_ylabels_align(h)
% IRF_PLOT_YLABELS_ALIGN left-align ylabels 
%
%   IRF_PLOT_YLABELS_ALIGN align ylabels of all subplots in current figure
%
%   IRF_PLOT_YLABELS_ALIGN(h) align ylabels of axis with handles h

narginchk(0,1)

if nargin==0 % get subplot handles and align all of them
    h=irf_plot_get_subplot_handles;
end

ylext=zeros(1,numel(h));
ylh=zeros(1,numel(h));
for jh=1:numel(h)
    ylh(jh)=get(h(jh),'ylabel');
    set(ylh(jh),'units','normalized'); % if other units, labels will move when zooming plot
    ext=get(ylh(jh),'extent');
    ylext(jh)=ext(1);
end
min_ylext=min(ylext);
for jh=1:numel(h)
    pos=get(ylh(jh),'position');
    pos(1)=pos(1)-(ylext(jh)-min_ylext);
    set(ylh(jh),'position',pos);
end
