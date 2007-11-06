function y = irf_gse2gei(x,flag_gse2gei)
% IRF_GSE2GEI  transform between GSE and GEI
%
% Y = IRF_GSE2GEI(X) 
%     Transforms X(GSE) -> Y(GEI)
%
% Y = IRF_GSE2GEI(X,-1) 
%     Transforms X(GEI) -> Y(GSE)
%
% GEI - Geocentric Equatorial Inertial
%   X - axis points from the Earth toward the first point of Aries 
%       (the position of Sun at the vernal equinox)
%   Z - axis is parallel to the rotation axis of the Earth and points
%   northward
%   Y - axis completes the right-handed orthogonal set Y = ZxX
%
% See also IRF_GSE2GSM
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(1,2,nargin))

if nargin==1, flag_gse2gei = 1; end

d = gei_gse_trans_matrix(x(1,1));
y = x;

if flag_gse2gei~=-1, d = d'; end    % GSE->GEI

for i=1:3
	y(:,i+1) = d(i,1).*x(:,2) + d(i,2).*x(:,3) + d(i,3).*x(:,4); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Help functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = tdt (jday)
% Correction factor to convert universal time (UT) to 
% terrestrial dynamic time (TDT)
res = (54.6 + 0.9 * (jday.integer - 2446066.50) / 365.) / 86400;
return

function y = fmod2p (x)
y = mod (x, 2*pi);
if y < 0, y = y + 2*pi; end
return

function y = fmod360 (x)
y = mod (x, 360);
if y < 0, y = y + 360; end
return

function res = kepler (e, m)
% Solves kepler - e*sin(kepler) = m
% e eccentricity
% m mean anomaly in degrees

int   i;
rm = m * pi/180; % in radians
k = rm;
for i = 1:10, k = rm + e * sin (k); end
res = k * 180 / pi;            % back to degrees */
return

function mat = rmatrix (angle, idx)
% input; rotation angle in degrees
% output: rotation  matrix
% index: 'X','Y','Z'

mat = zeros(3,3);

RAD = pi/180;
c = cos (RAD * angle);
s = sin (RAD * angle);

switch idx
	case 'X'                     % rotate about x axis
		mat([2 3],[2 3]) = [c, s; -s, c] ;
		mat(1,1) = 1;
    case 'Y'                     % rotate about y axis
		mat([1 3],[1 3]) = [c, -s; s, c] ;
		mat(2,2) = 1;
    case 'Z'                     % rotate about z axis
		mat([1 2],[1 2]) = [c, s; -s, c] ;
		mat(3,3) = 1;
	otherwise
		error('index must be X, Y or Z')
end
return

function rot_mat = rotate3 (angle1, index1, angle2, index2, angle3, index3)
% Returns compound rotation matrix
% angle in degrees, and axis, where index: 'X','Y','Z'

mat1 = rmatrix(angle1, index1);
mat2 = rmatrix(angle2, index2);
mat3 = rmatrix(angle3, index3);
rot_mat = ( mat1 * mat2 ) * mat3;
return

function res = ecliptic (jday)
t = (jday.integer - 2451545.) / 36525.; % rel to 2000.0
res = 23.439291 - 1.300417e-2 * t - 1.63889e-7 * t * t + 5.03611e-7 * t^3;
return

function c = crossn(a, b)
% Normalized vector product
c = cross(a, b);
x = sqrt(sum(c.*c));
if (x >= 1e-30), c = c / x; end
return

function geigse = gei_gse_trans_matrix(epoch)
e_days = fix( (epoch + 43200.0) /86400);
% J1970 = 2440587.50
jday.integer = 2440587 + e_days;
jday.fraction = epoch + 43200.0 - e_days*86400.0;

eqlipt = [0.0, -0.398, 0.917];
geigse = zeros(3,3);
geigse(:,1) = sun_vect(jday);
geigse(:,2) = crossn(eqlipt, geigse(:,1));
geigse(:,3) = crossn(geigse(:,1), geigse(:,2));

return

function geocentric = sun_vect(jday)
% Sun-Earth vector

eclip_v = earth(jday);
eclip_v = -eclip_v/sqrt(sum(eclip_v.*eclip_v)); % normalized earth -> sun

% ecliptic to equatorial
rot_mat = rotate3(-ecliptic (jday), 'X', 0, 'X', 0, 'X');

geocentric = eclip_v*rot_mat;

function vector2 = earth (jday)

lme = [ 279.69668, 36000.76892,.0003025, 0.; ...    % longitude
	358.47583, 35999.04975, -.000150, -.0000033; ...% mean anomaly
    .01675104, -.0000418, -.000000126, 0.];         % eccentricity

t = (tdt(jday) + (jday.integer - 2415020.0) + jday.fraction) / 36525.; % rel to 1900

t2 = t * t;
t3 = t2 * t;

elements = lme(:,1) + lme(:,2) * t + lme(:,3) * t2 + lme(:,4) * t3;

%  Additional perturbations, pg 82 Meeus
RAD = pi/180; AU = 1.4959965e8;
A = fmod2p (RAD * 153.23 + 22518.7541 * t); % venus
B = fmod2p (RAD * 216.57 + 45037.5082 * t); % venus
C = fmod2p (RAD * 312.69 + 32964.3557 * t); % jupiter
D = fmod2p (RAD * 350.74 + 445267.1142 * t -.00144 * t2); % moon
E = fmod2p (RAD * 231.19 + 20.2 * t);
H = fmod2p (RAD * 353.40 + 65928.7155 * t);


longitude = fmod360 (elements(1) + ...
	.00134 * cos (A) + ...
	.00154 * cos (B) + ...
	.00200 * cos (C) + ...
	.00179 * sin (D) + ...
	.00178 * sin (E));

Radius = 1.0000002 + ...
	.00000543 * sin (A) + ...
	.00002575 * sin (B) + ...
	.00001627 * sin (C) + ...
	.00003076 * cos (D) + ...
	.00000927 * sin (H);

mean_anomaly = fmod360 (elements(2));
eccentricity = elements(3);
mean_dist = Radius * AU;
%inclination = 0.;
%asc_node = 0.;

eccentric_anomaly = kepler (eccentricity, mean_anomaly);

vector1(1) = -mean_dist * (cos (RAD * eccentric_anomaly) - eccentricity);
vector1(2) = -mean_dist * ...
	(sqrt (1.- eccentricity * eccentricity) * sin (RAD * eccentric_anomaly));
vector1(3) = 0;

rot_mat = rotate3( mean_anomaly - longitude, 'Z', 0, 'X', 0, 'X');

vector2 = vector1*rot_mat;

return