function QV = thor_QV(R,rTHOR)
% Quality factor for bow shock crossing.
%   QV = THOR_QV(R,rTHOR);
%     R - bow shock nose standoff distance
%     rTHOR - coordinates of THOR
%
%   QV = cos(theta)^2
%     theta - angle between shock normal and GSE X, Sun direction

units = irf_units;

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
  
% Bowshock model
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
fx = @(y,R0) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+y.^2/0.04)+R0;

% First resample THORs location to the bowshock location
newX = fx(yTHOR,R0);
newY = fx(newX,R0);


x1 = newX-0.001; x2 = newY+0.001;
y1 = (fy(x1,R0)); y2 = (fy(x2,R0));
dotprod = [x2-x1]*1 + [y2-y1]*0;
normmat = repmat(sqrt((x2-x1).^2+(y2-y1).^2),1,1); 
alfa = acosd(dotprod./normmat);
theta = 90-alfa;

QV = real(cosd((theta)).^2);
%irf_plot({TSeries(Time,theta),TSeries(Time,real(QV)),TSeries(Time,imag(QV))})
if isa(rTHOR,'TSeries') || isa(R,'TSeries')
  QV = irf.ts_scalar(Time,QV); 
  QV.name = 'QV';
end