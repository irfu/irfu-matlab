function [hax,hcb] = plot_projection(varargin)
%MMS.PLOT_SKYMAP Plots projection.
%
% [ax,hcb] = MMS.PLOT_SKYMAP(TSeriesDesDist,'Opt1',OptVal1,...) - plot skymap.
% MMS.PLOT_SKYMAP(AX,...) - plot in axes AX.
%
% Options:
%   'tint' - plot data for time interval if tint.length = 2
%             or closest time if tint.length = 1
%   'vectors' - Nx2 cell array with 1x3 vector in first column and
%                textlabel in second column
%   'xyz' - 3x3 matrix with [x;y;z]. z is normal to the plotted plane and
%           x and y are made orthogonal to z and each other if they are 
%           not already. If you want to plot different planes you have to
%           rotate this matrix -> [y;z;x] -> [z;x;y]
%   'clim' - [cmin cmax], colorbar limits in logscale
%   'vlim' - vmax in km/s, zoom in to XLim = YLim = vlim*[-1 1]
%   'elevationlim' - in degrees [0 90], limiting elevation angle  
%                    above/below projection plane to include in projection
%
%   Example:     
%     mms.plot_projection(disDist1,'tint',tint+0.1*[-1 1],'xyz',[-1 0 0; 0 -1 0;hatB0],'vlim',1000);
%     mms.plot_projection(desDist1,'tint',tint,'xyz',[-1 0 0; 0 -1 0;hatB0],'vlim',40000);

irf.log('warning','Please verify that you think the projection is done properly!');

[ax,args,nargs] = axescheck(varargin{:});

dist = args{1}; args = args(2:end);
if isempty(dist); irf.log('warning','Empty input.'); return; end

doFlipX = 0;
doFlipY = 0;
have_vlabels = 0;
have_vectors = 0;
have_clim = 0;
have_vlim = 0;
limElevation = 20;
correctForBinSize = 0;

tId = 1:dist.length;
[~,xi] = hist([log10(10),log10(30e3)],32); energyTable = 10.^xi; % energy levels
energyEdges = 10.^linspace(log10(10),log10(30e3),33);

x = [1 0 0]; y = [0 1 0]; z = [0 0 1]; % default vectors

if nargs > 1, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))   
    case 'tint'
      l = 2;
      tint = args{2};
      if tint.length == 1
        tId = find( abs(dist.time-tint) == min(abs(dist.time-tint)) );        
      else
        [tId,~] = dist.time.tlim(tint); 
        if isempty(tId); irf.log('warning','No data for given time interval.'); return; end              
      end
    case 'vectors'
      l = 2;
      vectors = args{2};
      have_vectors = 1;
    case 'xyz'
      l = 2;
      coord_sys = args{2};
      x = coord_sys(1,:)/norm(coord_sys(1,:));
      y = coord_sys(2,:)/norm(coord_sys(2,:));
      z = coord_sys(3,:)/norm(coord_sys(3,:));
      y = cross(z,cross(y,z)); y = y/norm(y);
      x = cross(y,z); x = x/norm(x);
      if abs(acosd(y*(coord_sys(2,:)/norm(coord_sys(2,:)))'))>1; 
        irf.log('warning',['y (perp1) changed from [' num2str(coord_sys(2,:)/norm(coord_sys(2,:)),'% .2f') '] to [' num2str(y,'% .2f') '].']);
      end
      if abs(acosd(x*(coord_sys(1,:)/norm(coord_sys(1,:)))'))>1; 
        irf.log('warning',['x (perp2) changed from [' num2str(coord_sys(1,:)/norm(coord_sys(1,:)),'% .2f') '] to [' num2str(x,'% .2f') '].']);
      end
    case 'clim'
      l = 2;
      clim = args{2};
      have_clim = 1;
    case 'vlim'
      l = 2;
      vlim = args{2};
      have_vlim = 1;  
    case 'elevationlim'
      l = 2;
      limElevation = args{2}; 
    case 'vlabel'
      l = 2;
      have_vlabels = 1;
      vlabels = args{2};   
    case 'usebincorrection'      
      l = 2;
      correctForBinSize = args{2};
    case 'flipx'
      doFlipX = 1;
    case 'flipy'
      doFlipY = 1;  
  end
  args = args(l+1:end);  
  if isempty(args), break, end    
end

if have_vlabels
  vlabelx = vlabels{1};
  vlabely = vlabels{2};
  vlabelz = vlabels{3};
else
  vlabelx = ['v_{x=[' num2str(x,'% .2f') ']}'];
  vlabely = ['v_{y=[' num2str(y,'% .2f') ']}'];
  vlabelz = ['v_{z=[' num2str(z,'% .2f') ']}'];
end

% Check if it's electron or ions
if strfind(dist.name,'des')
  isDes = 1; 
elseif strfind(dist.name,'dis')
  isDes = 0;
else
  irf.log('warning','Can''t recognize if input is electron or ions. Assuming it''s electrons.');
  isDes = 1;
end

% Set up spherical coordinate system.
[~,polar] = hist([0 pi],16);
[~,azimuthal] = hist([0,2*pi],32);

[POL,AZ] = meshgrid(polar,azimuthal);
X = -sin(POL).*cos(AZ); % '-' because the data shows which direction the particles were coming from
Y = -sin(POL).*sin(AZ);
Z = -cos(POL);

% Transform into different coordinate system
xX = reshape(X,size(X,1)*size(X,2),1);
yY = reshape(Y,size(Y,1)*size(Y,2),1);
zZ = reshape(Z,size(Z,1)*size(Z,2),1);

newTmpX = [xX yY zZ]*x';
newTmpY = [xX yY zZ]*y';
newTmpZ = [xX yY zZ]*z';

newX = reshape(newTmpX,size(X,1),size(X,2));
newY = reshape(newTmpY,size(X,1),size(X,2));
newZ = reshape(newTmpZ,size(X,1),size(X,2));

elevationAngle = atan(newZ./sqrt(newX.^2+newY.^2));
planeAz = atan2(newY,newX)+pi;

% gets velocity in direction normal to 'z'-axis
geoFactorElev = cos(elevationAngle);
% geoFactorBinSize - detector bins in 'equator' plane are bigger and get a larger weight.
% I think this is not good for the implementation in this function
if correctForBinSize, geoFactorBinSize = sin(POL); 
else, geoFactorBinSize = 1; end

% New coordinate system
edgesAz = linspace(0,2*pi,33);
[nAz,binAz] = histc(planeAz,edgesAz);

% Collect data in new distribution function FF
FF = zeros(32,32); % azimuthal, energy
for iE = 1:32;
  for iAz = 1:32;
    % dist.data has dimensions nT x nE x nAz x nPol
    C = squeeze(irf.nanmean(dist.data(tId,iE,:,:),1));
    %C = squeeze(nansum(dist.data(tId,iE,:,:),1));
    C = C.*geoFactorElev.*geoFactorBinSize;    
    C(abs(elevationAngle)>limElevation*pi/180) = NaN;
    % surf(ax,newX,newY,newZ,C);xlabel(ax,'x');ylabel(ax,'y')
    C(planeAz<edgesAz(iAz)) = NaN;
    C(planeAz>edgesAz(iAz+1)) = NaN;   
    
    FF(iAz,iE) = irf.nanmean(irf.nanmean(C));
  end    
end
%%
% Plot projection plane
if isempty(ax), fig = figure; ax = axes; axis(ax,'square'); end

Units = irf_units; % Use IAU and CODATA values for fundamental constants.
if isDes, m = Units.me; else m = Units.mp; end
speedTable = sqrt(energyEdges*Units.e*2/m)*1e-3; % km/s
rE = speedTable;
%rE = energyTable;
plX = rE'*cos(edgesAz+pi);
plY = rE'*sin(edgesAz+pi);

if isDes; % make electron velocities 10^3 km/s
  hs = surf(ax,plX*1e-3,plY*1e-3,plY*0,log10(FF'*1e30)); % s+3*km-6
  vUnitStr= '(10^3 km/s)'; 
  if have_vlim, ax.YLim = vlim*[-1 1]*1e-3; ax.XLim = ax.YLim; end
  hcb = colorbar('peer',ax);
  hcb.YLabel.String = 'log_{10} f_e (s^3km^{-6})';
else % ion velocities km/s
  hs = surf(ax,plX,plY,plY*0,log10(FF'*1e30)); % s+3*km-6  
  vUnitStr= '(km/s)'; 
  if have_vlim, ax.YLim = vlim*[-1 1]; ax.XLim = ax.YLim; end
  hcb = colorbar('peer',ax);
  hcb.YLabel.String = 'log_{10} f_i (s^3km^{-6})';
end

ax.XLabel.String = [vlabelx ' ' vUnitStr];
ax.YLabel.String = [vlabely ' ' vUnitStr];
ax.ZLabel.String = [vlabelz ' ' vUnitStr];

axis(ax,'square')
axes(ax)
view([0 0 1])

if doFlipX
  cbPos = hcb.Position;
  axPos = ax.Position;
  view(ax,180,90)
  hcb.Position = cbPos;
  ax.Position = axPos;
end
ax.Box = 'on';

titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
ax.Title.String = titleString;   

shading(ax,'flat');
if have_clim, ax.CLim = clim; end

% Plot vectors
hold(ax,'on');

while have_vectors  
  if isempty(vectors), break, end  
  vecHat = vectors{1,1}/norm(vectors{1,1});
  vecHat = vecHat;
  vecTxt = vectors{1,2};
  scale = ax.YLim(2)/3;
  scaleQuiver = 1.5*scale;   
  quiver3(ax,0,0,0,(vecHat*x'),(vecHat*y'),abs(vecHat*z'),scaleQuiver,'linewidth',1,'linestyle','-','color','k')
  scaleTxt = 1.9*scale;
  text(double(scaleTxt*vecHat*x'),double(scaleTxt*vecHat*y'),1,vecTxt,'fontsize',10)
  vectors = vectors(2:end,:);  
  axis(ax,'square')
end
hold(ax,'off');

end
