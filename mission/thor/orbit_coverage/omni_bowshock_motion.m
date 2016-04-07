
%% define time interval
tintTemp = '2010-01-01T00:00:00/2010-01-11T00:00:00Z';
tint2013 = '2013-01-01T00:00:00/2013-12-31T23:59:00Z';

tintUTC = tint2013;
tint = irf.tint(tintUTC);

%% get omni pressure data
pbm_orig = irf_get_data(tint,'P,Bz,Ms','omni_min');

%% clean up data
pbm=pbm_orig;
pbm(isnan(pbm(:,2)),:)=[];
pbm(isnan(pbm(:,3)),:)=[];
pbm(isnan(pbm(:,4)),:)=[];
P = pbm(:,2);
Bz = pbm(:,3);
M = pbm(:,4);

%% calculate bow shock distance R
rzero=(10.22+1.29*tanh(0.184*(Bz+8.14))).*P.^(-1/6.6);
gamma=5/3;
M=5; % do not use OMNI because there are many NaNs
R =rzero.*(1+1.1*((gamma-1)*M.^2+2)./((gamma+1)*(M.^2-1)));

%% bin it along edges
dbin = 0.05; % RE
edges = 11:dbin:18;
[N,edges,binR] = histcounts(R,edges);

%% remove equal neighbours
isEqualNeighbour = (binR(1:end-1) == binR(2:end));
binR1 = binR;
binR1(isEqualNeighbour)=[];
disp([num2str(sum(isEqualNeighbour)) ' equal neigbours removed' ...
	' out of ' num2str(numel(binR)) ' points.']);

%% remove points passed by shock
binR2=binR1;
indPassedByShock = find(sign(binR1(3:end)-binR1(2:end-1)) == ...
	sign(binR1(2:end-1)-binR1(1:end-2)))+1;
binR2(indPassedByShock)=[];
disp([num2str(numel(indPassedByShock)) ' points passed by shock removed' ...
	' out of ' num2str(numel(binR1)) ' points.']);

%% count crossings
nCrossings = zeros(1,size(edges,2)-1);
for ii = 2:numel(binR2)
	indShocksCrossed = binR2(ii-1)+1 : binR2(ii);
	nCrossings(indShocksCrossed) = nCrossings(indShocksCrossed) + 1;
end

%% statistics of shock crossings [#/day]
totalTimeDays = (tint.stop.tts-tint.start.tts)/(24*3600);
probabilityCrossings = nCrossings / totalTimeDays;
distanceRE = (edges(2:end)+edges(1:end-1))/2;

% probability matrix, first column distance, second column #/day
% probabilityCrossings = [distanceRE' nCrossings' / totalTimeDays];

%% statistics of shock distance 
countsR=histcounts(R,edges);


%% plot the probability of shock crossing
figure;
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

hb=bar(distanceRE,probabilityCrossings,1);
grid on;
xlabel('Distance [R_E]');
ylabel('Probability of shock crossing [#/day]');

%% add the probability of shock didstance 

hold(hb.Parent,'on');
hb2=bar(hb.Parent,distanceRE,countsR/100,1,'facealpha',0.5,'facecolor','green');
xax = get(hb.Parent,'xlim');
yax = get(hb.Parent,'ylim');
ht=text(xax(1)*1.01,yax(1)*1.01,'Probability of BS distance',...
	'color','green','rotation',90);

%% mark bow shock box requirement

hold(hb.Parent,'on');
Rmin = 13;
Rmax = 15;
yax = get(hb.Parent,'ylim');
plot(hb.Parent,[Rmin Rmin],yax,'r','linewidth',2);
plot(hb.Parent,[Rmax Rmax],yax,'r','linewidth',2);
title(hb.Parent,tintUTC);
text(Rmin*1.01,yax(2)*0.95,'Bow shock ROI','color','red')


%% ESTIMATE fraction of quasi-parallel shock crossings and quasi perp
%% get omni pressure data
bbbm = irf_get_data(tint,'Bx,By,Bz,Ms','omni_min');
%% clean up data
om=bbbm;
for i=2:size(om,2)
	om(isnan(om(:,i)),:)=[];
end
Bx = om(:,2);
By = om(:,3);
Bz = om(:,4);
Ms = om(:,5);

%% Calculate angle
angleDeg = atand(sqrt(By.^2+Bz.^2)./abs(Bx));

%% Plot angle probability
figure;
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

h=histogram(angleDeg,'normalization','cdf');
set(h.Parent,'ylim',[0 1]);
title(h.Parent,'Shock angle probability at the BS nose')
text(0,0.95,tintUTC)
grid on;
xlabel('Shock angle [deg]');
ylabel('Cumulative density function');

