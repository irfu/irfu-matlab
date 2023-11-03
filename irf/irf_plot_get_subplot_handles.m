function h=irf_plot_get_subplot_handles(figurehandle)
%IRF_PLOT_GET_SUBPLOT_HANDLES return subplot handles of current figure

if nargin==0, figurehandle=gcf;end

ud=get(figurehandle,'userdata');
if isfield(ud,'subplot_handles')
  h=ud.subplot_handles;
  h=h(ishandle(h)); % remove h values that are not handles
  if numel(h)~=numel(ud.subplot_handles) % in case some h values had to be removed update ud
    ud.subplot_handles=h;
    set(figurehandle,'userdata',ud);
  end
else

  htmp=findobj(figurehandle,'type','axes','-not','tag','Colorbar');
  number_of_subplots=numel(htmp);

  if number_of_subplots>=1 % if there are subplots
    axy=zeros(1,number_of_subplots);
    axx=axy;
    for ih=1:length(htmp)
      ax=get(htmp(ih),'position');
      axy(ih)=ax(2);
      axx(ih)=ax(1);
    end
    ind_ax=find(axx<0.2); % find
    hax=htmp(ind_ax);
    [~,ind]=sort(axy(ind_ax));
    ind=fliplr(ind);
    h=hax(ind);
  else
    h=[];
  end
end


