function irf_legend(axis_handle,labels,position,varargin)
% irf_legend(axis_handle,labels,position,text_property,text_value,...)
%
% labels - cell array with strings
% position - in normalized units 
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

cluster_colors=[[0 0 0];[1 0 0];[0 1 0];[0 0 1]];
colord=get(axis_handle, 'ColorOrder');

for i=1:length(labels),
  ht=text(position(1),position(2),labels{i},'parent',axis_handle,'units','normalized','verticalalignment','baseline','fontweight','bold','fontsize',12);
  set(ht,'color',colord(i,:));
  for j=1:size(varargin,2)/2
      textprop=varargin{2*j-1};
      textvalue=varargin{2*j};
      if strcmpi(textprop,'color') && strcmp(textvalue,'cluster') && i<=4,
          set(ht,'color',cluster_colors(i,:));  
      else
          set(ht,varargin{2*j-1},varargin{2*j});
      end
  end
  ext=get(ht,'extent'); position(1)=position(1)+ext(3)*1.4;
end