% Four ways to visualize data using color coding.
% Questions? Contact konrad.steinvall@irfu.se

%% Load example data to plot
tint=irf.tint('2017-07-26T07:01:20.00Z','2017-07-26T07:02:10.00Z');

sc=3;
c_eval('Bgse=mms.get_data(''B_gse_fgm_brst_l2'',tint+[-1,1],?);',sc);
c_eval('Egse=mms.get_data(''E_gse_edp_brst_l2'',tint+[-1,1],?);',sc);
c_eval('Bscm=mms.get_data(''B_gse_scm_brst_l2'',tint+[-1,1],?);',sc);
Efac = irf_convert_fac(Egse,Bgse,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm,Bgse,[1 0 0]);

Efacf=Efac.filt(50,0,[],5);
Bscmfacf=Bscmfac.filt(50,0,[],5);

%% #1: 2D histogram
xdata=Bscmfacf.x.data; %Bperp2
ydata=Efacf.x.data; %Eperp1

Nx = 200; %number of grid points in the x-dimension
Ny = 200; %number of grid points in the y-dimension

xedges = linspace(min(xdata),max(xdata),Nx); %min() and max() can be replaced
%by desired max and min values to improve plot. Necessary when there are
%large outliers.
yedges = linspace(min(ydata),max(ydata),Ny);

figure;
histpl=histogram2(xdata,ydata,xedges,yedges,'DisplayStyle','Tile');
histpl.EdgeColor='none';
maxcounts = max(max(histpl.BinCounts));
ax=gca;
cmap = colormap(parula(maxcounts));
cmap(1,:)=[0.8,0.8,0.8]; %<- Set 1-count level to grey
cmap(2,:)=[0.8,0.8,0.8]; %<- Set 1-count level to grey
colormap(ax,cmap);
cbar=colorbar;
ylabel(cbar,'Counts');
xlabel('X-variable');
ylabel('Y-variable');

%% #2: Color-coded scatterplot

xdata=Bgse.x.data; %Bperp2
ydata=Bgse.y.data; %Eperp1
cdata = Bgse.time-Bgse.time(1); % Color-code seconds since first data point.
%Can also be any other quantity with the same nr of datapoints as xdata and
%ydata.
cdata2 = Bgse.abs.data; %Magnitude of B.

dotsize=50;
figure;
fig=gcf;
fig.Position(3:4)=[1000,400];
subplot(1,2,1)
scatter(xdata,ydata,dotsize,cdata,'filled');
cbar=colorbar;
ylabel(cbar,'Seconds since t_0')
xlabel('B_x (nT)');
ylabel('B_y (nT)');

subplot(1,2,2)
scatter(xdata,ydata,dotsize,cdata2,'filled');
cbar=colorbar;
ylabel(cbar,'|B| (nT)')
xlabel('B_x (nT)');
ylabel('B_y (nT)');


%% #3: Color-coded lineplot (Nicer, but more complicated than #2
% Figure cannot be saved in vector format ('painters').

shortTint1=irf.tint('2017-07-26T07:01:46.2726Z','2017-07-26T07:01:46.2841Z');
shortTint3=shortTint1+[0.002,-0.002];

figure;
fig=gcf;
fig.Position(3:4)=[500,400];
subplot(1,1,1);
ax=gca;
testpl=plot(Bscmfacf.tlim(shortTint3).x.data,Bscmfacf.tlim(shortTint3).y.data,'linewidth',4);
xlabel('\deltaB_{\perp1} (nT)' ,'interpreter','tex');
ylabel('\deltaB_{\perp2} (nT)' ,'interpreter','tex');

% Setup the colormap
datapoints = length(Bscmfacf.tlim(shortTint3).time);
cdata = [uint8(parula(datapoints)*255) uint8(ones(datapoints,1))].';
% Apply the colormap to the data
drawnow
set(testpl.Edge, 'ColorBinding','interpolated', 'ColorData',cdata)

cbar2=colorbar(ax);
%Put colorbar above the plot:
cbar2.Location='northoutside';
cbar2.Ticks=[];
cbarlabel = ylabel(cbar2,'Time $\longrightarrow$','interpreter','latex','Fontsize',16);


%% #4: Color-coded TSeries (slightly messy)
% This code plots a normal TSeries (inputTS), but the color of the TSeries can
% depends on the values of a second TSeries (colorTS) sharing the same
% timeline, i.e. resampled to the same timeline.

% Example 1:
% This is perhaps not the best example, because Bscm changes a lot, so the
% plot looks discontinuous. Another example is shown below.

h=irf_plot(1,'newfigure');
fig=gcf;
fig.Position(3:4)=[1200,300]; % Stretch fig as it would look in a normal
% plot with several panels.

% We want to color Bscm based on Epar in a short interval.
Tshort = irf.tint('2017-07-26T07:01:46.613495Z','2017-07-26T07:01:46.772801Z');

colorTS = Efacf.z.resample(Bscmfacf.tlim(Tshort));

% We need to specify the "bins" of the color. In this case, E|| goes from
% roughly -100 to 100, so we set
% Could for example use min(colorTS), max(colorTS), but outliers may mess
% that up.
vals=-100:20:100;
[newTS,colors]=irf_plot_linecolor_2TS(Bscmfacf.x.tlim(Tshort),colorTS,vals);
hold(h(1),'on');
for ii=1:length(newTS)
  irf_plot(h(1),newTS{ii},'o','color',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',5,'linewidth',1)
  if ii<length(newTS)
    % This part is just to fix the colorbar labels which can probably
    % be done in a better way.
    charvec{ii+1,1}=num2str(vals(ii));
  else
    charvec{ii+1}=[];
  end
end
for itest=3:2:length(charvec)-2

  charvec{itest}='';
end
ax=h(1);
colormap(ax,colors)
cbar=colorbar(ax);
cbar.Ticks = [0:1/size(colors,1):1];
cbar.TickLabels = charvec;
ylabel('B_{scm,x} (nT)','interpreter','tex');
ylabel(cbar,'E_{||} (mV/m)','fontsize',22);

%Set background to grey to make it easier to see the data
h(1).Color=[0.75,0.75,0.75];


%% #4 Example 2: Color Vix by |B|


c_eval('Vi = mms.get_data(''Vi_gse_fpi_brst_l2'',tint+[-10,10],?);',sc);
c_eval('Bgse_long=mms.get_data(''B_gse_fgm_brst_l2'',tint+[-10,10],?);',sc);


h=irf_plot(1,'newfigure');
fig=gcf;
fig.Position(3:4)=[1200,300]; % Stretch fig as it would look in a normal
% plot with several panels.

colorTS = Bgse_long.abs.resample(Vi);
% Change Vals according to values of colorTS
% Could for example use min(colorTS), max(colorTS), but outliers may mess
% that up.
vals=10:1:25;
[newTS,colors]=irf_plot_linecolor_2TS(Vi,colorTS,vals);
hold(h(1),'on');
for ii=1:length(newTS)
  irf_plot(h(1),newTS{ii},'o','color',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',6,'linewidth',1)
  if ii<length(newTS)
    charvec{ii+1,1}=num2str(vals(ii));
  else
    charvec{ii+1}=[];
  end
end
for itest=3:2:length(charvec)-2

  charvec{itest}='';
end
ax=h(1);
colormap(ax,colors)
cbar=colorbar(ax);
cbar.Ticks = [0:1/size(colors,1):1];
cbar.TickLabels = charvec;
ylabel('Vi{x} (km/s)','interpreter','tex');
ylabel(cbar,'|B| (nT)','fontsize',22);

%Set background to grey to make it easier to see the data
h(1).Color=[0.75,0.75,0.75];