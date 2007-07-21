function h=irf_plot_get_subplot_handles
%IRF_PLOT_GET_SUBPLOT_HANDLES return subplot handles of current figure

  htmp=findobj(gcf,'type','axes','-not','tag','Colorbar');
  hmax=1;
  for ih=1:length(htmp),
    ax=get(htmp(ih),'position');
    axy(ih)=ax(2);axx(ih)=ax(1);
  end
  ind_ax=find(axx<0.2); % find 
  hax=htmp(ind_ax);
  [xx,ind]=sort(axy(ind_ax));
  ind=fliplr(ind);
  h=hax(ind);

