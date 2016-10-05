function [ax,hcb] = plot_skymap(varargin)
% MMS.PLOT_SKYMAP Plots skymap.
%
%  [ax,hcb] = MMS.PLOT_SKYMAP(TSeriesDesDist,'Opt1',OptVal1,...);
%     ax - handle to axes
%     hcb - handle to colorbar
%  MMS.PLOT_SKYMAP(AX,...); - plot in axes AX.
%
%  Options:
%    'tint' - plot data for time interval if tint.length = 2
%             or closest time if tint.length = 1
%    'energy' - plot energy closest to given energy
%    'energylevel' - plot given energylevel
%    'vectors' - Nx2 cell array with 1x3 vector in first column and
%                textlabel in second column, eg. 
%                vectors = {Bhat,'B';Ehat,'E'}
%    'flat' - plot a flat skymap (ie. not a sphere)
%    'log' - plot log10 scale
%    'energytable' - energytable from v1 data
%    'phi' - phi data froom V1 data, vary from time to time, ax.XLim = [-10 370]; 
%    'phib' & 'polarb' [no theta] - solid angle boundary, only for 'flat'; 

[ax,args,nargs] = axescheck(varargin{:});

dist = args{1}; args = args(2:end);
if isempty(dist); irf.log('warning','Empty input.'); return; end

plotLog = 0; fString = 'f (s^3km^{-6})';
plotSphere = 1;
plotb = 0;          % flag for plotting boundary; 2016-03-11, wyli
have_vectors = 0;
tId = 1:dist.length;
eId = 1:32;
[~,xi] = hist([log10(10),log10(30e3)],32); energyTable = 10.^xi;
% Set up spherical coordinate system.
r = 1; % radius of sphere
phi_edges = linspace(0,2*pi,size(dist.data,3)+1);  % azimuthal angle bin edges, default
theta_edges = linspace(0,pi,size(dist.data,4)+1); % polar angle bin edges, default

if nargs > 1, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))
    case 'energytables'
        l = 2;
        energyTable = args{2};       
    case 'energy'
      l = 2;
      energy = args{2};      
      eId = find(abs(energyTable-energy)==min(abs(energyTable-energy)));
    case 'energylevel'
      l = 2;
      eId = args{2};
    case 'phi_edges'            % phi edges from V1/V2 data
      l = 2;
      phi_edges = args{2};
    case 'phib'            % phi boundary for picking partial distribution
      l = 2;
      phib = args{2};
      plotb = 1;            % flag --> 1
    case 'polarb'            % theta boundary for picking partial distribution
      l = 2;
      polarb = args{2};      
    case 'tint'
      l = 2;
      tint = args{2};
      if tint.length == 1
        tId = find(abs(dist.time-tint)==min(abs(dist.time-tint)));        
      else
        [tId,~] = dist.time.tlim(tint);       
        if isempty(tId); irf.log('warning','No data for given time interval.'); return; end
      end
    case 'vectors'
      l = 2;
      vectors = args{2};
      have_vectors = 1;
    case 'flat'
      plotSphere = 0;
    case 'log'
      plotLog = 1;
      fString = 'log_{10} f (s^3km^{-6})';
  end
  args = args(l+1:end);
  if isempty(args), break, end  
end

[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles were coming from
Y = -r*sin(THETA).*sin(PHI);
Z = -r*cos(THETA);
% dist.data has dimensions nT x nE x nAz x nPol
C = squeeze(nanmean(nanmean(dist.data(tId,eId,:,:),2),1))';


% Plot skymap
if isempty(ax), fig = figure; ax = axes; end
if plotLog, fC = log10(C*1e30); % s+3*km-6 
else fC = (C*1e30); end % s+3*km-6; 
  
if plotSphere
  hs = surf(ax,X,Y,Z,fC); % 
  axis(ax,'square')
  axis(ax,'equal')
  ax.XLabel.String = 'X';
  ax.YLabel.String = 'Y';
  ax.ZLabel.String = 'Z';
  %titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(electronEnergyLevels(k))]};  
  shading(ax,'flat');
else % plot flat map
  [~,theta] = hist(theta_edges,16);
  [~,phi] = hist(phi_edges,32);   
  hs = surf(ax,PHI*180/pi,THETA*180/pi,THETA*0,[fC(:,17:32) fC(:,1:16)]);  
  ax.XLabel.String = 'Azimuthal angle (deg)';
  ax.YLabel.String = 'Polar angle (deg)';
  shading(ax,'flat');
  view(ax,[0 0 1])  
  ax.YLim = [0 180];
  ax.XLim = [-10 370];
  ax.Box = 'on';
end
hcb = colorbar('peer',ax);
hcb.YLabel.String = fString;
titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'],['Energy = ' num2str(energyTable(eId),'%.0f') ,' eV']};
ax.Title.String = titleString;   

% Plot vectors
hold(ax,'on');

while have_vectors
  if plotSphere
    vecHat = vectors{1,1}/norm(vectors{1,1});
    vecTxt = vectors{1,2};

    scale = 1.5;   
    quiver3(ax,-scale*vecHat(1),-scale*vecHat(2),-scale*vecHat(3),vecHat(1),vecHat(2),vecHat(3),2*scale,'linewidth',2)
    scale = 1.7;
    axes(ax)
    text(double(scale*vecHat(1)),double(scale*vecHat(2)),double(scale*vecHat(3)),vecTxt,'fontsize',14)    
  else % plot flat skymap
    vecHat = vectors{1,1}/norm(vectors{1,1});
    vecTxt = vectors{1,2};
    [azim,elev,r] = cart2sph(vecHat(1),vecHat(2),vecHat(3));
    if azim<0, azim = azim + 2*pi; end
    if azim>2*pi, azim = azim - 2*pi; end
    
    plot(ax,azim*180/pi,elev*180/pi+90,'o','linewidth',2,'markersize',12,'color',[1 0 0])
    plot(ax,azim*180/pi,elev*180/pi+90,'o','linewidth',0.5,'markersize',2,'color',[1 0 0],'markerfacecolor',[0 0 0])
    axes(ax); text(double(azim*180/pi),double(elev*180/pi+90),['   ' vecTxt],'fontsize',14,'HorizontalAlignment','left', 'color', [1 1 1])
    
    [azim,elev,r] = cart2sph(-vecHat(1),-vecHat(2),-vecHat(3)); 
    if azim<0, azim = azim + 2*pi; end
    if azim>2*pi, azim = azim - 2*pi; end
    plot3(ax,azim*180/pi,elev*180/pi+90,0,'o','linewidth',2,'markersize',12,'color',[1 0 0])
    plot3(ax,azim*180/pi,elev*180/pi+90,0,'x','linewidth',2,'markersize',12,'color',[1 0 0])       
    axes(ax); text(double(azim*180/pi),double(elev*180/pi+90),['   ' vecTxt],'fontsize',14,'HorizontalAlignment','left', 'color', [1 1 1])        
  end  
  vectors = vectors(2:end,:);
  if isempty(vectors), break, end  
end
    
  if plotb
    if not(plotSphere)            
    phib = [phib, phib(1)];
    polarb = [polarb, polarb(1)];        
    plot(ax, phib, polarb, 'linewidth',2, 'color',[1 1 1])
    end
  end  
hold(ax,'off');

end