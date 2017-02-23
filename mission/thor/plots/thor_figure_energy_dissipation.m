

xb=28; % breakpoint between lines
close
hca=axes; hold(hca,'on')
set(gcf,'position',[560 639 560 309]);

fontsize=16;
set(gcf,'defaultAxesFontSize',fontsize);
set(gcf,'defaultTextFontSize',fontsize);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
% lines
x1 = 5:xb; 
x2 = xb:40;
y0 = 100;
k1 = -1; 
k2 = -2;
y1 = y0+k1*x1;
y2 = y1(end)*(1-k2/k1)+k2/k1*y0+k2*x2;
x =[x1 x2];
y =[y1 y2];
plot(hca,x,y,'linewidth',2)
ylabel(hca,'log(Energy density)','fontsize',fontsize)
xlabel(hca,'log(1/length scale)','fontsize',fontsize)
set(gca,'xtick',[],'ytick',[],...
    'xlim',5*[-1 1]+x([1 end]),'ylim',5*[-2 1]+sort(y([1 end])))

% arrow, energy flow
arx0=10;
ary0=84;
arlength=10; 
arrow([arx0 ary0],[arx0+arlength ary0+k1*(arlength)],'width',3,'tipangle',35)
artext=text(arx0-3,ary0-8,'Energy flow','rotation',-0);

% dissipation region
plot([xb xb],[42 85],'--k')
text(32,82,'THOR','rotation',-0,'fontsize',18,'color',[0 0 0]);
text(32,74,{'dissipation','range'},'rotation',-0);

% scales
y_scales=45;
if 1    
    %arrow([xb-11 y_scales],[xb-15 y_scales],'width',2,'tipangle',35,'baseangle',60)
    %arrow([xb+11 y_scales],[xb+15 y_scales],'width',2,'tipangle',35,'baseangle',60)
    arrow([xb-1 y_scales],[xb-13 y_scales],'width',2,'tipangle',35,'baseangle',60)
    arrow([xb+1 y_scales],[xb+13 y_scales],'width',2,'tipangle',35,'baseangle',60)
    text(xb-1,y_scales+3,{'fluid scale'},'rotation',-0,'horizontalalignment','right');
    th=text(xb+1,y_scales+3,'kinetic scale','rotation',-0);
else
    arrow([xb-3 y_scales],[xb y_scales],'width',2,'tipangle',35,'baseangle',60)
    arrow([xb+3 y_scales],[xb y_scales],'width',2,'tipangle',35,'baseangle',60)
    text(xb-13,y_scales,{'larger scale'},'rotation',-0);
    text(xb+4,y_scales,'kinetic scale','rotation',-0);
end
box(gca,'on')

