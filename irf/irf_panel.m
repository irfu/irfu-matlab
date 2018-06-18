function hout=irf_panel(handles,tag)
%IRF_PANEL  set current subplot axis
%
%   HCA = IRF_PANEL(subplotnumber);
%   HCA = IRF_PANEL(tag);
%   HCA = IRF_PANEL(handles,tag);
%   HCA = IRF_PANEL(handles,subplotnumber);
%   handles - handles of figures subplots
%   tag - unique string (for example current date) identifying subplot
%   subplotnumber  - number of subplot

flag_tag_defined=0;
if nargin==1 && ischar(handles)
    tag=handles;
    if isempty(get(0,'CurrentFigure')) % there is no figure open
        irf_figure(1);              % create new figure with one panel
    end
    ud=get(gcf,'userdata');
    if isfield(ud,'subplot_handles') && any(ishandle(ud.subplot_handles))
        handles=ud.subplot_handles;
    else
		handles=gca;
		ud.subplot_handles=handles;
		set(gcf,'userdata',ud);
    end
    flag_tag_defined=1;
end
h=handles;
parent=get(h(1),'parent');
ud=get(parent,'userdata');
subplot_handles=ud.subplot_handles;
    
if nargin==1 && ~flag_tag_defined 
    h=handles;
    parent=get(h(1),'parent');
    ud=get(parent,'userdata');
    if isfield(ud,'current_panel')
        hout=ud.current_panel;
    else
        hout=handles(1);
        ud.current_panel=hout;
        set(parent,'userdata',ud);
    end
else
    if ischar(tag) % check the tag
        hca=findobj(h,'tag',tag);
        if numel(hca)>0 % has found subplot with tag
            hout=hca(1);
            irf.log('warning',['--SUBPLOT-- <' tag '> (Using existing panel)'])
            hnumber=find(hout==subplot_handles);
            parent=get(hout,'parent');
            ud=get(parent,'userdata');
            ud.current_panel=hnumber;
            set(parent,'userdata',ud);
        else % go to next subplot and add tag to subplot
            parent=get(h(1),'parent');
            ud=get(parent,'userdata');
            if isfield(ud,'current_panel')
                current_panel=ud.current_panel+1;
                if current_panel > numel(h)
                    current_panel=numel(h);
                end
            else
                current_panel=1;
            end
            hout=h(current_panel);
            set(hout,'tag',tag);
            irf.log('warning',['--SUBPLOT-- <' tag '> (New panel)'])
            ud.current_panel=current_panel;
            set(parent,'userdata',ud); 
        end
    elseif isnumeric(tag) % set subplot number
        hout=handles(tag);
        parent=get(hout,'parent');
        ud=get(parent,'userdata');
        ud.current_plot=tag;
        set(parent,'userdata',ud);
				if nargout==0
					axes(hout);
				end
    end    
end
