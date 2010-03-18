function irf_legend(axis_handle,labels,position,varargin)
% irf_legend(labels,position,axis_handle,text_property,text_value,...)
%
% labels - cell array with strings
% position - in normalized units 
%
% irf_legend({'B_X','B_Y','B_Z','B'},[0 

colord=get(axis_handle, 'ColorOrder');

for i=1:length(labels),
  ht=text(position(1),position(2),labels{i},'parent',axis_handle,'units','normalized','verticalalignment','baseline','fontweight','bold','fontsize',12);
  set(ht,'color',colord(i,:));
  for j=1:size(varargin,2)/2
      set(ht,varargin{2*j-1},varargin{2*j});
  end
  ext=get(ht,'extent'); position(1)=position(1)+ext(3)*1.2;
end