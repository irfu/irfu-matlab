function irf_legend(axis_handle,labels,position,varargin)
% irf_legend(axis_handle,labels,position,text_property,text_value,...)
%
% labels - cell array with strings
% position - in normalized units 
%
% Examples:
% irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02, 0.1]);
%
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[0 0 0];[1 0 0]]);
% irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02, 0.1]);


colord=get(axis_handle, 'ColorOrder');

for i=1:length(labels),
  ht=text(position(1),position(2),labels{i},'parent',axis_handle,'units','normalized','verticalalignment','baseline','fontweight','bold','fontsize',12);
  set(ht,'color',colord(i,:));
  for j=1:size(varargin,2)/2
      set(ht,varargin{2*j-1},varargin{2*j});
  end
  ext=get(ht,'extent'); position(1)=position(1)+ext(3)*1.4;
end