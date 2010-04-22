% IRFNOTES File with different common examples how to use irf routines
 
edit irfnotes; return

%% Initializing some figure
for ha=1, % define size to have best agreement with eps file
set(0,'defaultLineLineWidth', 1.5);
fn=figure(61);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
end
for ha=1, % set subplots
% specifying position
h(1)=axes('position',[0.65 0.78 0.2 0.2]); % [x y dx dy]
% having all in standard form
n_subplots=8;i_subplot=1;
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% 

end
%% Add information to figures
for ha=1, % text and legends
ht=irf_pl_info([mfilename '  ' datestr(now)]);set(ht,'interpreter','none');
end
for ha=1, % labels a),b)...
numb={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
for ip=1:2,
  axes(h(ip));
  ht=irf_pl_info(numb{ip},gca,[0.01,1]);
  set(ht,'fontsize',10,'verticalalignment','top');
end
end
%% Second axis 
for ha=1, % example
hl1 = line(x1,y1,'Color','r');
ax1 = gca;
set(ax1,'XColor','r','YColor','r')

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
end
