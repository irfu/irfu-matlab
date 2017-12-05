clf
h = plot3([-1 1],[0 0],[0 0],'k'); %HFA-X
h.LineWidth=3;
hold on
h = plot3([0 0],[-1 1],[0 0],'r'); %HFA-Y
h.LineWidth=3;
h = plot3([0 0],[0 0],[-1 1],'g'); %HFA-z
h.LineWidth=3;

sc = 1.5; off = 0.5;
h =plot3(sc*[-1 1]+off,sc*[-1 1]+off,sc*[-1 1]+off,'b'); %Boom
h.LineWidth=3;

h =quiver3(0,0,0,1,-1,1); % Sun
h.LineWidth=3;
legend('HFA-X','HFA-Y','HFA-Z','boom','Sun')

%%

clf

h = plot3([-1 1],[-1 1],[-1 1],'k') %HFA-X
h.LineWidth=3;
hold on
h = plot3([-1 1],[1 -1],[-1 1],'r') %HFA-Y
h.LineWidth=3;
h = plot3([-1 1],[-1 1],[1 -1],'g') %HFA-Z
h.LineWidth=3;

h =quiver3(0,0,0,1,0,0);

sc = 1.5;
h =plot3([0 0],sc*[-1 3],[0 0],'b');
h.LineWidth=3;

xlabel('X - towards sun')
ylabel('Y - towards sc')
zlabel('Z')

l = 3;
set(gca,'XLim',l*[-1 1],'YLim',l*[-1 1],'ZLim',l*[-1 1])