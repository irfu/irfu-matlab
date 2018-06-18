function irf_colorbar_fit_label_height(hcb)
% IRF_COLORBAR_FIT_LABEL_HEIGHT fit label within the colorbar hcb
hcbl = get(hcb,'ylabel');
hy=hcbl;
colorbar_label_fontsize=get(hy,'fontsize');
units=get(hy,'units');
set(hy,'units','normalized');
temp=get(hy,'Extent');
colorbarlabelheight = temp(4);
while colorbarlabelheight>1.1
  colorbar_label_fontsize=colorbar_label_fontsize*0.95;
  %set(hy,'fontsize',colorbar_label_fontsize,'position',labelposition);
  set(hy,'fontsize',colorbar_label_fontsize);
  temp=get(hy,'Extent');
  colorbarlabelheight=temp(4);
end

end
