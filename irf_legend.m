function ht=irf_legend(axis_handle,labels,position,varargin)
% irf_legend(axis_handle,labels,position,text_property,text_value,...)
%
% labels - cell array with strings
% position - in normalized units (default, if x position is epoch assumes data units)
%
% default for normalized units is "smart" alignment,e.g.'left' in left side 
%  and 'right' in righ side of plot
%
% nonstandard text property values: 
%  'color'='cluster' - if 4 labels given write them in cluster colors (black,red,green,blue)
%
% Examples:
% irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02, 0.1]);
%
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[0 0 0];[1 0 0]]);
% irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02, 0.1]);
%
% irf_legend(gca,{'C1','C2','C3','C4'},[0.02, 0.9],'color','cluster')

%hgca = gca;

if axis_handle == 0, % add to the whole figure if inhandle = 0
    axis_handle = axes('Units','normalized', 'Position',[0 0 1 1], 'Visible','off', ...
        'Tag','BackgroundAxes', 'HitTest','off');
end

unit_format='normalized';
colord=get(axis_handle, 'ColorOrder');
cluster_colors=[[0 0 0];[1 0 0];[0 0.5 0];[0 0 1]];

% use smart alignment in upper part of panel vertical alignment='top', 
% when vertical position above 1 then again 'bottom', in bottom part use
% vertical alignment 'bottom'. Similarly for horizontal alignment

if position(1)>1e8 % assume that position specifued in axis units
    unit_format='data';
    ud = get(gcf,'userdata');
    if isfield(ud,'t_start_epoch')
        position(1)=position(1)- double(ud.t_start_epoch);
    end
end

if position(1)<0, 
    value_horizontal_alignment='right';
elseif position(1) < 0.5 
    value_horizontal_alignment='left';
elseif position(1) <= 1 
    value_horizontal_alignment='right';
else    
    value_horizontal_alignment='left';
end

if position(2)<0, 
    value_vertical_alignment='top';
elseif position(2) < 0.5 
    value_vertical_alignment='baseline';
elseif position(2) <= 1 
    value_vertical_alignment='top';
else    
    value_vertical_alignment='baseline';
end


if ischar(labels), % Try to get variable labels from string (space separates).
    lab=labels;clear labels;
    labels{1}=lab;
    colord=[0 0 0];
end

if strcmpi(value_horizontal_alignment,'left'), 
    for i=1:length(labels), % start with first label first
        ht(i)=text(position(1),position(2),labels{i},'parent',axis_handle,'units',unit_format,'fontweight','demi','fontsize',12);
        set(ht(i),'color',colord(i,:));
        set(ht(i),'verticalalignment',value_vertical_alignment);
        set(ht(i),'horizontalalignment',value_horizontal_alignment);
        for j=1:size(varargin,2)/2
            textprop=varargin{2*j-1};
            textvalue=varargin{2*j};
            if strcmpi(textprop,'color') && strcmp(textvalue,'cluster') && i<=4,
                set(ht(i),'color',cluster_colors(i,:));
            else
                set(ht(i),varargin{2*j-1},varargin{2*j});
            end
        end
        ext=get(ht(i),'extent'); position(1)=position(1)+ext(3)*1.4;
    end
else
    for i=length(labels):-1:1, % start with last label first
        ht(i)=text(position(1),position(2),labels{i},'parent',axis_handle,'units','normalized','fontweight','demi','fontsize',12);
        set(ht(i),'color',colord(i,:));
        set(ht(i),'verticalalignment',value_vertical_alignment);
        set(ht(i),'horizontalalignment',value_horizontal_alignment);
        for j=1:size(varargin,2)/2
            textprop=varargin{2*j-1};
            textvalue=varargin{2*j};
            if strcmpi(textprop,'color') && strcmp(textvalue,'cluster') && i<=4,
                set(ht(i),'color',cluster_colors(i,:));
            else
                set(ht(i),varargin{2*j-1},varargin{2*j});
            end
        end
        ext=get(ht(i),'extent'); position(1)=position(1)-ext(3)*1.4; %% !!! minus sign because we start with last label and go left
    end    
end

