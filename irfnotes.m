% IRFNOTES File with different common examples how to use irf routines
%
% use code folding options to fast find your necessary examples

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
%% Reading files
for ha=1, % formatted file reading
%File contents are time intervals in format "T1 T2 Comments":
%2008-03-03T22:50:00 2008-03-03T23:30:00
%2008-03-10T22:10:00 2008-03-10T22:45:00 !
%2008-03-13T07:40:00 2008-03-13T09:40:00 ? shock?

%file reading 
[t1,t2,tint_comments]=textread('Events_reconnection.txt','%s%s%[^\n]');
for j=1:size(t1,1),
  tint(j,1)=iso2epoch(t1{j});tint(j,2)=iso2epoch(t2{j});
end
clear t1 t2 j;
end
%% Cluster data reading
for ha=1, % using c_get_batch
    c_get_batch(toepoch([2002 03 04 10 00 00]),30*60,'sp','/home/yuri/caa-data/20020304')
% if time intervals to download are in matrix tint 
for j=1:size(tint,1),
 c_get_batch(tint(j,1),tint(j,2)-tint(j,1),'sp',['./' epoch2iso(tint(j,1),1) '-' epoch2iso(tint(j,2),1)]);
end
clear j;
end

%% Other
