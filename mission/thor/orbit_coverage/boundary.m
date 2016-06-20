function [x,y] = boundary(R0,boundary)
% BOUNDARY Gives bowshock or magnetopause coordinates.
%   [x,y] = boundary(R0,boundary);
%       R0 - Standoff distance in RE
%       boundary - 'mp' gives magnetopause, also default
%                  'bs' gives bowshock
%       x - GSE x coordinate in RE
%       y - GSE y coordinate in RE

if isempty(boundary); boundary = 'mp'; end

% Magnetopause
% Reference: Shue et al 1998
% Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
% Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
% Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
% Default values: Dp=2nPa, Bz=0nT
Bz = 0; Dp = 2;

theta=0:0.1:pi;
alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
r=R0*(2./(1+cos(theta))).^alpha;
xMP=r.*cos(theta);
yMP=r.*sin(theta);
yMP(abs(xMP)>100)=[];
xMP(abs(xMP)>100)=[];

% Bowshock
xBS=R0:-0.5:-100;
yBS=sqrt(0.04*(xBS-R0).^2-45.3*(xBS-R0)); % original F/G model adds rstandoff^2=645

switch lower(boundary)
    case 'bs'; x=xBS; y=yBS;
    case 'mp'; x=xMP; y=yMP;
end

x = [x x];
y = [y -y];
[y,isort] = sort(y);
x = x(isort);

