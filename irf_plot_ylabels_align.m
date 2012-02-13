function irf_plot_ylabels_align(h)
% IRF_PLOT_YLABELS_ALIGN left-align ylabels 
%
%   IRF_PLOT_YLABELS_ALIGN align ylabels of all subplots in current figure
%
%   IRF_PLOT_YLABELS_ALIGN(h) align ylabels of axis with handles h

error(nargchk(0,1,nargin))

if nargin==0, % get subplot handles and align all of them
    h=irf_plot_get_subplot_handles;
end

ylext=zeros(1,numel(h));
for jh=1:numel(h)
    ylh=get(h(jh),'ylabel');
    un=get(ylh,'units');
    set(ylh,'units','normalized');
    ext=get(ylh,'extent');
    ylext(jh)=ext(1);
    set(ylh,'units',un);
end
min_ylext=min(ylext);
for jh=1:numel(h)
    ylh=get(h(jh),'ylabel');
    un=get(ylh,'units');
    set(ylh,'units','normalized');
    pos=get(ylh,'position');
    pos(1)=pos(1)-(ylext(jh)-min_ylext);
    set(ylh,'position',pos);
    set(ylh,'units',un);
end
