function varargout = thor_QB(R,rTHOR,B,varargin)
% Quality factor for bow shock crossing.
%   QB = THOR_QB(R,rTHOR,B);
%     R - bow shock nose standoff distance
%     rTHOR - coordinates of THOR
%     B - interplanetary magnetic field
%
%     QV = cos(theta)^2
%     theta - angle between shock normal and magnetic field

units = irf_units;

doParShock = 1;
if nargin>3 && strcmp(varargin{1},'perp')
  doParShock = 0;
end
  
if isa(R,'TSeries'); 
  Time = R.time;
  R0 = R.data; 
  if strcmp(R.units,'km') || R0(1,1) > 50; 
    R0 = R0/units.RE*1e3; % km->RE
  end  
else 
  R0 = R;
end
if isa(rTHOR,'TSeries'); 
  Time = rTHOR.time;
  xTHOR = rTHOR.x.data;
  yTHOR = rTHOR.y.data;
  if strcmp(rTHOR.units,'km') || rTHOR.data(end,1) > 100; 
    xTHOR = xTHOR/units.RE*1e3; % km->RE
    yTHOR = yTHOR/units.RE*1e3; % km->RE
  end
else
  xTHOR = rTHOR(:,1);
  yTHOR = rTHOR(:,2);
end
if isa(B,'TSeries'); 
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
if 0 % Normalize magnetic field
  Bnorm = sqrt(Bx.^2+By.^2+0*Bz.^2);
  Bx = Bx./Bnorm;
  By = By./Bnorm;
  Bz = Bz./Bnorm;
end

  
% Bowshock model
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
fx = @(y,R0) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+y.^2/0.04)+R0;

% First resample THORs location to the bowshock location
newX = fx(yTHOR,R0);
newY = fy(newX,R0);

x1 = newX-0.0001; x2 = newX;
y1 = (fy(x1,R0)); y2 = (fy(x2,R0));
y1(yTHOR<0) = -y1(yTHOR<0);
y2(yTHOR<0) = -y2(yTHOR<0);

dx12 = [x2-x1];
dy12 = [y2-y1];

dnorm = sqrt(dx12.^2+dy12.^2);

x3 = x1 + dx12*cosd(90)-dy12*sind(90);
y3 = y1 + dx12*sind(90)+dy12*cosd(90);

dnorm = sqrt(dx12.^2+dy12.^2);
dx13 = [x3-x1];
dy13 = [y3-y1];
dnorm13 = sqrt(dx13.^2+dy13.^2);
dx13 = dx13./dnorm13;
dy13 = dy13./dnorm13;
%aa=acosd(dx23.*dx12+dy23.*dy12);

normdir = irf.ts_vec_xyz(Time,[dx13,dy13,dy13*0]); tsTheta.name = 'Bowshock surface';

dotprod = dx13.*Bx + dy13.*By + 0.*Bz;
%normmat = repmat(sqrt(dx13.^2+dx13.^2),1,1); 
theta = acosd(dotprod);%./normmat);%dotprod./normmat);
%alfaXY = acosd(dx13.*Bx + dy13.*By + 0.*Bz);%dotprod./normmat);
%alfaZ = atand(Bz./sqrt(Bx.^2+By.^2));

%alfa=alfaXY;
tsTheta = irf.ts_scalar(Time,theta); tsTheta.name = 'Magnetic field normal angle';
if doParShock;
  %theta = alfa;
  QB = real(cosd((theta)).^2);
else
  %theta = alfa;
  QB = real(sind((theta)).^2);
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
elseif nargout == 2;
  %tsTheta = irf.ts_scalar(Time,theta); tsTheta.name = 'Magnetic field normal angle';
  varargout = {QB,tsTheta};
elseif nargout == 3;
  %tsTheta = irf.ts_scalar(Time,theta); tsTheta.name = 'Magnetic field normal angle';
  varargout = {QB,tsTheta,normdir};  
end