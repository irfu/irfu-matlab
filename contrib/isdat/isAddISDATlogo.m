function isAddISDATlogo()
%isAddISDATlogo - Puts an ISDAT logo onto a plot.
islogo=imread('isPlotStamp.png','BackgroundColor', [1 1 1]);
axes('position',[0.9 0 0.1 0.05]);
image(islogo), set(gca,'Visible','off');
