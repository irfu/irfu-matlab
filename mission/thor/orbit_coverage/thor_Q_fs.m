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

xRef = 15;
yRef = 15;
%vecRef = [xRef yRef]./norm([xRef yRef]);
vecTHOR_ = ([xTHOR yTHOR]);
vecTHOR = irf_norm(vecTHOR_);
vecRef_ = [xTHOR-xRef yTHOR-yRef];
vecRef = irf_norm(vecRef_);


dotProdArg = vecRef(1)*vecTHOR(:,1)+vecRef(2)*vecTHOR(:,2);
dotProd = dotProdArg;
thetaTHOR = acosd(dotProd);

tsTheta      = irf.ts_scalar(Time,abs(thetaTHOR)); 
tsTheta.name = 'Angle from [5 -5]RE';
tsTheta.units = 'deg';

Q = cosd(thetaTHOR).^2;
% Different Q, just distance from Y=0
%Q = 0.5-0.5*tanh((abs(yTHOR)-10)/10);
Q = exp(-abs(yTHOR).^2/20^2);
tsQ = irf.ts_scalar(Time,Q); 
tsQ.name = 'Q fs';



if nargout == 1
  varargout = {tsQ};
elseif nargout == 2
  varargout = {tsQ,tsTheta};
end