function iKSR = thor_in_ksrs(rTHOR,B,Dp,M,BSNX)
% Finds the times THOR spends in different key science regions (KSR's)
%
% [iKSR] = thor_isCrossing(rTHOR,B,Dp,M);
%   Input: B - solar wind magnetic field (nT)
%          Dp - solar wind dynamic pressure (nP/m2)
%          M - solar wind Mach number
%   Output:
%   iKSR = 1 - magnetosphere, inside magnetopause
%          2 - magnetosheath, outside magnetopause but inside bowshock
%          3 - foreshock, outside bowshock but inside the
%              region limited by the tangent of IMF B to the bowshock
%          4 - pristine solar wind, outside bowshock and foreshock
%          5 - magnetopause crossing
%          6 - bowshock crossing
%



% Magnetopause model, Shue 1998
theta=0:0.1:pi;
f_mpR0 = @(Bz,Dp) (10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6);
f_alpha = @(Bz,Dp) (0.58-0.007*Bz)*(1+0.024*log(Dp));
f_mpR = @(Bz,Dp,theta) f_mpR0(Bs,Dp)*(2./(1+cos(theta))).^f_alpha(Bz,Dp);

% Bowshock model
gamma = 5/3;
f_bsR0 = @(Bz,Dp,M) f_mpR(Bz,Dp).*(1+1.1*((gamma-1).*M.^2+2)./((gamma+1)*(M.^2-1)));
f_bsR = @(x,R0) sqrt(0.04*(x-f_bsR0(Bz,Dp,M)).^2-45.3*(x-f_bsR0(Bz,Dp,M))); % original F/G model adds rstandoff^2=645
f_bsX = @(r,R0) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+r.^2/0.04)+f_bsR0(Bz,Dp,M);


%x=r.*cos(theta);
%y=r.*sin(theta);
%ii=find(abs(x)>100);
%x(ii)=[];
%y(ii)=[];


%% Check if THOR is inside bowshock or not
iKSR = irf.ts_scalar(rTHOR.time,nan(rTHOR.length,1));
useConstantBforForeshock = 1;
Bspiral = [-1,1,0];

[outsideMP,insideMP,crossingMP] = magnetopause(rTHOR,B.z.data,Dp.data);
[outsideBS,insideBS,crossingBS] = bowshock(rTHOR,B.z.data,Dp.data,M.data,BSNX.data);
if useConstantBforForeshock
  [outsideFS,insideFS,crossingFS] = foreshock(rTHOR,Bspiral,BSNX.data);
else
  [outsideFS,insideFS,crossingFS] = foreshock(rTHOR,B.data,BSNX.data);
end


%plot(rTHOR(outsideMP).x.data,rTHOR(outsideMP).y.data,'.r',rTHOR(insideMP).x.data,rTHOR(insideMP).y.data,'.b',rTHOR(crossingMP).x.data,rTHOR(crossingMP).y.data,'.g')
%plot(rTHOR(outsideBS).x.data,rTHOR(outsideBS).y.data,'.r',rTHOR(insideBS).x.data,rTHOR(insideBS).y.data,'.b',rTHOR(crossingBS).x.data,rTHOR(crossingBS).y.data,'.g')
%plot(rTHOR(isOutside).x.data,rTHOR(isOutside).y.data,'or',rTHOR(isInside).x.data,rTHOR(isInside).y.data,'.b')

allInd = 1:rTHOR.length;
isMSP = setdiff(allInd,outsideMP);
isMSH = intersect(outsideMP,insideBS);
isFS = intersect(outsideBS,insideFS);
isPSW = intersect(outsideBS,outsideFS);

%plot(rTHOR(isMSP).x.data,rTHOR(isMSP).y.data,'.r')
%plot(rTHOR(isMSH).x.data,rTHOR(isMSH).y.data,'.r')
iKSR.data(isMSP) = ones(numel(isMSP),1);
iKSR.data(isMSH) = repmat(2,numel(isMSH),1);
iKSR.data(isFS) = repmat(3,numel(isFS),1);
iKSR.data(isPSW) = repmat(4,numel(isPSW),1);
iKSR.data(crossingMP) = repmat(5,numel(crossingMP),1);
iKSR.data(crossingBS) = repmat(6,numel(crossingBS),1);

iKSR.userData.indexKSR.i_1 = 'magnetosphere, inside magnetopause';
iKSR.userData.indexKSR.i_2 = 'magnetosheath, outside magnetopause but inside bowshock';
iKSR.userData.indexKSR.i_3 = 'foreshock, outside bowshock but inside the region limited by the tangent of IMF B to the bowshock';
iKSR.userData.indexKSR.i_4 = 'pristine solar wind, outside bowshock and foreshock';
iKSR.userData.indexKSR.i_5 = 'magnetopause crossing';
iKSR.userData.indexKSR.i_6 = 'bowshock crossing';

if useConstantBforForeshock
  iKSR.userData.indexKSR.i_3 = sprintf('foreshock, outside bowshock but inside the region limited by the tangent of IMF B/|B| = [%.1f,%.1f,%.1f]; to the bowshock',Bspiral/norm(Bspiral));
end

return


end
% functions
function [isOutside,isInside,isCrossing] = magnetopause(rTHOR,B,Dp)
units = irf_units;
xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
yTHOR = rTHOR.y.data/units.RE*1e3;
zTHOR = rTHOR.z.data/units.RE*1e3;
r_THOR = sqrt(xTHOR.^2+yTHOR.^2+zTHOR.^2);

% B = Bz
% Magnetopause model
f_mpR0 = @(B,Dp) (10.22+1.29*tanh(0.184*(B+8.14))).*Dp.^(-1/6.6); % magnetopause standoff distance
mpR0 = f_mpR0(B,Dp);
%mpR0 = repmat(12,size(mpR0));
f_alpha = @(B,Dp) (0.58-0.007*B).*(1+0.024*log(Dp));
f_mpR = @(B,Dp,theta) f_mpR0(B,Dp).*(2./(1+cos(theta))).^f_alpha(B,Dp);


%alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
thetaTHOR = acos(xTHOR./sqrt(xTHOR.^2+yTHOR.^2+zTHOR.^2));
mpR = f_mpR(B,Dp,thetaTHOR);
%     x=mpR.*cos(theta);
%     y=mpR.*sin(theta);

allInd = 1:rTHOR.length;
isInside = find(xTHOR<mpR0 & abs(mpR)>abs(r_THOR));
isOutside = tocolumn(setdiff(allInd,isInside));

isInboundCrossing = isInside(find(diff(isInside)>1));
isOutboundCrossing = isOutside(find(diff(isOutside)>1));
isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
end

function [isOutside,isInside,isCrossing] = bowshock(rTHOR,B,Dp,M,BSNX)
units = irf_units;
xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
yTHOR = rTHOR.y.data/units.RE*1e3;
zTHOR = rTHOR.z.data/units.RE*1e3;
r_THOR = sqrt(yTHOR.^2+zTHOR.^2);

% Magnetopause model
f_mpR0 = @(Bz,Dp) (10.22+1.29*tanh(0.184*(Bz+8.14))).*Dp.^(-1/6.6); % magnetopause standoff distance
mpR0 = f_mpR0(B,Dp);

% Bowshock model
gamma = 5/3;
f_bsR0 = @(Bz,Dp,M) f_mpR0(Bz,Dp).*(1+1.1*((gamma-1).*M.^2+2)./((gamma+1)*(M.^2-1))); % bowshock standoff distance
if 0
  bsR0 = f_bsR0(B,Dp,M);
  f_bsR = @(x,Bz,Dp,M) sqrt(0.04*(x-f_bsR0(Bz,Dp,M)).^2-45.3*(x-f_bsR0(Bz,Dp,M))); % original F/G model adds rstandoff^2=645
  %f_bsX = @(r,Bz,Dp,M) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+r.^2/0.04)+f_bsR0(Bz,Dp,M);
  bsR0 = f_bsR0(B,Dp,M);
  rBS = f_bsR(xTHOR,B,Dp,M);
else
  bsR0 = BSNX;
  f_bsR = @(x,bsR0) sqrt(0.04*(x-bsR0).^2-45.3*(x-bsR0)); % original F/G model adds rstandoff^2=645
  rBS = (f_bsR(xTHOR,bsR0));
end
%rBS(isimag(rBS)) = ;
%tsrBS = irf.ts_scalar(mpR0.time,rBS);

allInd = 1:rTHOR.length;
isInside = find(xTHOR<bsR0 & abs(rBS)>abs(r_THOR));
isOutside = tocolumn(setdiff(allInd,isInside));

isInboundCrossing = isInside(find(diff(isInside)>1));
isOutboundCrossing = isOutside(find(diff(isOutside)>1));
isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
end

function [isOutside,isInside,isCrossing] = foreshock(rTHOR,B,BSNX)
units = irf_units;
xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
yTHOR = rTHOR.y.data/units.RE*1e3;
zTHOR = rTHOR.z.data/units.RE*1e3;
r_THOR = sqrt(yTHOR.^2+zTHOR.^2);

B = B./norm(B);
B0 = B;
%B0 = [-1,1,0]; % parker spiral
% Bowshock model
BAng = atan2d(B(:,2),B(:,1));
BAng = atan2d(B0(:,2),B0(:,1));
theta = linspace(0,180,30); % 'polar' angle, from x_GSE
azimuthal = linspace(0,360,50);

bsR0 = BSNX;
bsR0mean = mean(bsR0);
scaleR = bsR0/bsR0mean;
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
% use the mean value and then just rescale with BSNX
x = bsR0mean*cosd(theta);
y = fy(x,bsR0mean); % original F/G model adds rstandoff^2=645

r = sqrt(x.^2+y.^2);
theta = atan2d(y,x)*pi/180;



% Transform to 3D
azim = pi/180*azimuthal; % Azimuthal angle [0 360]
elev = pi/2-theta; % Angle with respect to the ZY-plame, pi/2 is towards Sun
[AZ,ELEV] = meshgrid(azim,elev);

% Bowshock surface
[Y,Z,X] = sph2cart(AZ,ELEV,repmat(r,numel(azim),1)');
% Forshock normal
[nX, nY, nZ] = surfnorm(X,Y,Z);
% Find point where magnetic field is tangent to surface
tangentNormal = cross(B0/norm(B0),[0 0 1]); % 45 deg to parker spiral, in xy-plane
tangNormalAngle = acosd(-(nX*tangentNormal(1) + nY*tangentNormal(2) + nZ*tangentNormal(3)));
indTangent = find(tangNormalAngle==min(min(tangNormalAngle)),1,'first');
% rescale X to bowshock distance
scaleX = 1;%bsR0/bsR0mean;
tX = X(indTangent)*scaleX; % RE
tY = Y(indTangent)+scaleX*0;
tZ = Z(indTangent)+scaleX*0;
tAng = tangNormalAngle(indTangent)*180/pi;

rTHOR_ = rTHOR-[tX tY tZ]*units.RE*1e-3; rTHOR_ = rTHOR_.data;
% 0 deg towards Sun (x_GSE)
angleXGSE = atan2d(rTHOR_(:,2),rTHOR_(:,1));
FS = rTHOR;

indNotFS = find(angleXGSE>BAng-180);
indNotFS2 = find(angleXGSE<-BAng-180);
allInd = 1:rTHOR.length;
isOutside = unique([indNotFS;indNotFS2]);
isInside = tocolumn(setdiff(allInd,isOutside));


isInboundCrossing = isInside(find(diff(isInside)>1));
isOutboundCrossing = isOutside(find(diff(isOutside)>1));
isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
%FS.data(indNotFS) = NaN;
%FS.data(indNotFS2) = NaN;
%  FS.data(indBS) = NaN;

1;
end
