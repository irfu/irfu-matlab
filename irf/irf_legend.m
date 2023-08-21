function ht=irf_legend(varargin)
% IRF_LEGEND add legend or text to axis
%
% IRF_LEGEND(AX,..) add to the specified axes
%
% IRF_LEGEND(labels,position,text_property,text_value,...)
%       add labels to the given position
%
% labels - Cell array with strings. Labels are aligned horizontally if
%           labels is a row vector and vertically if labels is a column
%           vector.
% position - in normalized units (default, if x position is epoch assumes data units)
%
% default for normalized units is "smart" alignment,e.g.'left' in left side
%  and 'right' in righ side of plot
%
% nonstandard text property values:
%  'color'='cluster' - if 4 labels given write them in cluster colors (black,red,green,blue)
%  'color'='mms' - MMS colors
%
% Examples:
% irf_legend({'B_X','B_Y','B_Z','B'},[0.02, 0.1]);
%
% irf_legend({'B_X';'B_Y';'B_Z';'B'},[1.02, 0.75]); % vertical alignment
%
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[0 0 0];[1 0 0]]);
% irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02, 0.1]);
%
% irf_legend({'C1','C2','C3','C4'},[0.02, 0.9],'color','cluster')

[axis_handle,args,nargs] = axescheck(varargin{:});

if nargs == 0 % show only help
  help irf_legend;
  return
end

if isempty(axis_handle)
  if any(ishandle(args{1})) % first argument is axis handles
    if args{1}==0 % add to the whole figure
      axis_handle = axes('Units','normalized', 'Position',[0 0 1 1], 'Visible','off', ...
        'Tag','BackgroundAxes', 'HitTest','off');
      uistack(axis_handle,'bottom') % move axis to background
    else
      axis_handle=args{1};
    end
    args=args(2:end);
    nargs=nargs-1;
  else % no axis handles
    if isnumeric(args{1}) %
      if args{1}==0 % add to the whole figure
        axis_handle = axes('Units','normalized', 'Position',[0 0 1 1], 'Visible','off', ...
          'Tag','BackgroundAxes', 'HitTest','off');
        uistack(axis_handle,'bottom') % move axis to background
        args=args(2:end);
        nargs=nargs-1;
      else
        axis_handle=gca;
      end
    else
      axis_handle=gca;
    end
  end
end

if nargs<2
  error('IRFU_MATLAB:irf_legend:InvalidNumberOfInputs','Incorrect number of input arguments')
else
  labels=args{1};
  position=args{2};
  pvpairs=args(3:end);
end

unit_format='normalized';
switch class(axis_handle)
  case 'matlab.graphics.axis.Axes'
    colord = get(axis_handle, 'ColorOrder');
  case 'matlab.graphics.chart.primitive.Line'
    for ii = 1:numel(axis_handle)
      colord(ii,:) = axis_handle(ii).Color;
      axis_handle(ii) = axis_handle(ii).Parent;
    end
    axis_handle = unique(axis_handle);
end
cluster_colors=[[0 0 0];[1 0 0];[0 0.5 0];[0 0 1];[1 0 1];[1 1 0];[0 1 1]];
mms_colors=[[0 0 0];[213,94,0];[0,158,115];[86,180,233]]/255;


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

if position(1)<0
  value_horizontal_alignment='right';
elseif position(1) < 0.5
  value_horizontal_alignment='left';
elseif position(1) <= 1
  value_horizontal_alignment='right';
else
  value_horizontal_alignment='left';
end

if position(2)<0
  value_vertical_alignment='top';
elseif position(2) < 0.5
  value_vertical_alignment='baseline';
elseif position(2) <= 1
  value_vertical_alignment='top';
else
  value_vertical_alignment='baseline';
end


if ischar(labels) % Try to get variable labels from string (space separates).
  lab=labels;clear labels;
  labels{1}=lab;
  colord=[0 0 0];
end

% If labels are given in a row vector, they should be plotted horizontally
% If labels are given in a column vector, they should be plotted vertically
labsize = size(labels);
if labsize(1)>1 && labsize(2)==1 % column vector
  value_alignment = 'vertical';
else % row vector
  value_alignment = 'horizontal';
end

% print labels in the order 1,2,...,N if left-aligned horizontal or top-aligned vertical
if (strcmp(value_alignment,'horizontal') && strcmpi(value_horizontal_alignment,'left')) || ...
    (strcmp(value_alignment,'vertical') && strcmpi(value_vertical_alignment,'top'))
  label_order=1:length(labels);
else % otherwise print labels in the order N,N-1,...,1
  label_order=length(labels):-1:1;
end
ht=gobjects(1,length(labels)); % allocate handles
tmp_ref_pos=position(1); % reference position in x
tmp_ref_ext_y=0; % if vertical add this to ht(i).Postion(2)
% loop through labels
for i=label_order % start with first label first
  ht(i)=text(position(1),position(2),labels{i},'parent',axis_handle,'units',unit_format,'fontweight','normal');
  set(ht(i),'color',colord(i,:));
  set(ht(i),'verticalalignment',value_vertical_alignment);
  set(ht(i),'horizontalalignment',value_horizontal_alignment);
  % loop through options ('color','FontSize','Interpreter',...)
  for j=1:size(pvpairs,2)/2
    textprop=pvpairs{2*j-1};
    textvalue=pvpairs{2*j};
    if strcmpi(textprop,'verticalalignment')
      value_vertical_alignment=textvalue; % value has been reset manually by input parameter
    end
    if strcmpi(textprop,'horizontalalignment')
      value_horizontal_alignment=textvalue; % value has been reset manually by input parameter
    end
    if strcmpi(textprop,'color') && strcmp(textvalue,'cluster') && i<=4
      set(ht(i),'color',cluster_colors(i,:));
    elseif strcmpi(textprop,'color') && strcmp(textvalue,'mms') && i<=4
      set(ht(i),'color',mms_colors(i,:));
    else
      set(ht(i),textprop,textvalue);
    end
  end
  % Get position and extent of label just printed
  txt_ext=get(ht(i),'extent');
  txt_pos=get(ht(i),'position');
  %
  if strcmp(value_alignment,'horizontal') % if horizontal labels
    if strcmpi(value_horizontal_alignment,'left')
      txt_pos(1)=txt_pos(1)-(txt_ext(1)-tmp_ref_pos); % how much to shift wrt to previous label
      set(ht(i),'position',txt_pos);
      txt_ext=get(ht(i),'extent');
      tmp_ref_pos=txt_ext(1)+txt_ext(3)+txt_ext(3)/max(1,numel(labels{i})); % the new reference position
    elseif strcmpi(value_horizontal_alignment,'right')
      txt_pos(1)=txt_pos(1)-(txt_ext(1)+txt_ext(3)-tmp_ref_pos);
      set(ht(i),'position',txt_pos);
      txt_ext=get(ht(i),'extent');
      tmp_ref_pos=txt_ext(1)-txt_ext(3)/max(1,numel(labels{i}));
    end
  elseif strcmp(value_alignment,'vertical') % if vertical labels
    ht(i).Position(2) = ht(i).Position(2)+tmp_ref_ext_y; % add to y-position
    % update what should be added to the y-position next time (can be < 0)
    if strcmp(value_vertical_alignment,'baseline')
      tmp_ref_ext_y = tmp_ref_ext_y+ht(i).Extent(4);
    elseif strcmp(value_vertical_alignment,'top')
      tmp_ref_ext_y = tmp_ref_ext_y-ht(i).Extent(4);
    end
  end

end

if nargout==0, clear ht; end

