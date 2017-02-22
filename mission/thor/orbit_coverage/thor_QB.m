function varargout = thor_QB(R,rTHOR,B,varargin)
% THOR_QB Quality factor for bow shock crossings
%
% [QB,tsTheta,normdir] = THOR_QB(R,rTHOR,B);
%   R     - bow shock nose standoff distance
%   rTHOR - coordinates of THOR
%   B     - interplanetary magnetic field
%
%   QB    - quality factor, default QB = cos(theta)^2 corresponding to
%           quasi-parallel shock
%   theta - angle between shock normal and magnetic field
%
% [QB,tsTheta,normdir] = THOR_QB(R,rTHOR,B,'perp');
%  estimate quality factor for quasi-perpendicular shock, where
%  QB=sin(theta)^2

units = irf_units;

doParShock = true;
if nargin>3 && strcmp(varargin{1},'perp')
  doParShock = false;
end
  
if isa(R,'TSeries')
  Time = R.time;
  R0 = R.data; 
  if strcmp(R.units,'km') || R0(1,1) > 50
    R0 = R0/units.RE*1e3; % km->RE
  end  
else 
  R0 = R;
end
if isa(rTHOR,'TSeries')
  Time = rTHOR.time;
%  xTHOR = rTHOR.x.data; %% NOT USED
  yTHOR = rTHOR.y.data;
  zTHOR = rTHOR.z.data;
  if strcmp(rTHOR.units,'km') || rTHOR.data(end,1) > 100
%    xTHOR = xTHOR/units.RE*1e3; % km->RE %% NOT USED
    yTHOR = yTHOR/units.RE*1e3; % km->RE
    zTHOR = zTHOR/units.RE*1e3; % km->RE
  end
else
%  xTHOR = rTHOR(:,1);  %% NOT USED
  yTHOR = rTHOR(:,2);
  zTHOR = rTHOR(:,3);
end
if isa(B,'TSeries')
  Time = B.time;
  B = B/B.abs;
  Bx = B.x.data;
  By = B.y.data;
  Bz = B.z.data;
else
  Bx = B(:,1);
  By = B(:,2);
  Bz = B(:,3);
end
if 1 % Normalize magnetic field
  Bnorm = sqrt(Bx.^2+By.^2+Bz.^2);
  Bx = Bx./Bnorm;
  By = By./Bnorm;
  Bz = Bz./Bnorm;
end

% THOR position in spherical coordinates
thetaTHOR = atand(zTHOR./yTHOR);
rTHOR = sqrt(yTHOR.^2+zTHOR.^2);

% Bowshock model
% Spherical coordinates, r,z (z=x)
fr = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0));
fx = @(r,R0) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+r.^2/0.04)+R0;

% First resample THORs location to the bowshock location
newX = fx(rTHOR,R0);
newR = fr(newX,R0);
newY = newR.*cosd(thetaTHOR);
newZ = newR.*sind(thetaTHOR);

% Using three points: first one is S/C, second one is a shift in the x direction
% and third one is a shift in the azimuthal direction.
x1 = newX; y1 = newY; z1 = newZ; r1 = newR;

x2 = x1 - 0.0001; r2 = fr(x2,R0); 
y2 = r2.*cosd(thetaTHOR); 
z2 = r2.*sind(thetaTHOR);

theta3 = thetaTHOR - 0.001; x3=x1;
y3 = r1.*cosd(theta3); 
z3 = r1.*sind(theta3);

y1(yTHOR<0) = -abs(y1(yTHOR<0));
y2(yTHOR<0) = -abs(y2(yTHOR<0));
y3(yTHOR<0) = -abs(y3(yTHOR<0));
z1(zTHOR<0) = -abs(z1(zTHOR<0));
z2(zTHOR<0) = -abs(z2(zTHOR<0));
z3(zTHOR<0) = -abs(z3(zTHOR<0));

% Creating two vectors (on the BS surface) from these three points.
d21 = [x1-x2, y1-y2, z1-z2]; d23 = [x3-x2, y3-y2, z3-z2];

% Cross product gives a vector perpendicular to the surface.
cr213 =  cross(d21,d23);

% Normalizing to unit vector.
normat = sqrt(cr213(:,1).^2+cr213(:,2).^2+cr213(:,3).^2);
vecx = cr213(:,1)./normat;
vecy = cr213(:,2)./normat;
vecz = cr213(:,3)./normat;

% Normal vector to the shock's surface at the spacecraft location.
nvec = [vecx, vecy, vecz];
normdir = irf.ts_vec_xyz(Time,[nvec(:,1),nvec(:,2),nvec(:,3)]);

% Dot product between normalized B field vector and normal vector.
dotprod = nvec(:,1).*Bx + nvec(:,2).*By + nvec(:,3).*Bz;
theta = acosd(dotprod);

tsTheta      = irf.ts_scalar(Time,theta); 
tsTheta.name = 'Magnetic field normal angle';

if doParShock
  QB = real(cosd((theta)).^2).*(theta<45); % QB=0 for angles corresponding to quasi-perp shock
else
  QB = real(sind((theta)).^2).*(theta>=45);% QB=0 for angles corresponding to quasi-par shock
end

if isa(rTHOR,'TSeries') || isa(R,'TSeries')
  QB = irf.ts_scalar(Time,QB); 
  if doParShock
    QB.name = 'QB par';
  else
    QB.name = 'QB perp';
  end
end

if nargout == 1
  varargout = {QB};
elseif nargout == 2
  varargout = {QB,tsTheta};
elseif nargout == 3
  varargout = {QB,tsTheta,normdir};  
end
