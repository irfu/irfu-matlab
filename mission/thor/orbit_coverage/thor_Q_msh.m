function varargout = thor_Q_msh(rTHOR)
% THOR_QB Quality factor for bow shock crossings
%
% [Q,tsTheta] = THOR_Q_MSH(R,rTHOR);
%   rTHOR - coordinates of THOR
%
%   Q    - quality factor, Q = cos(theta)^2
%   theta - angle between Sun-Earth line and spacecraft

units = irf_units;

Time = rTHOR.time;
xTHOR = rTHOR.x.data; 
yTHOR = rTHOR.y.data;
zTHOR = rTHOR.z.data;
if strcmp(rTHOR.units,'km') || rTHOR.data(end,1) > 100
  xTHOR = xTHOR/units.RE*1e3; % km->RE
  yTHOR = yTHOR/units.RE*1e3; % km->RE
  zTHOR = zTHOR/units.RE*1e3; % km->RE
end

% THOR in spherical coordinates
thetaTHOR = atan2d(sqrt(zTHOR.^2+yTHOR.^2),xTHOR);
tsTheta      = irf.ts_scalar(Time,abs(thetaTHOR)); 
tsTheta.name = 'Angle from bowshock nose';
tsTheta.units = 'deg';

% Quality factor, midnight is Q = 0; terminators are 0.7071
thetaQ = abs(thetaTHOR); thetaQ(thetaQ>90) = 90;
Q = cosd(thetaQ*0.5).^2;
tsQ = irf.ts_scalar(Time,Q); 
tsQ.name = 'Q msh';


if nargout == 1
  varargout = {tsQ};
elseif nargout == 2
  varargout = {tsQ,tsTheta};
end