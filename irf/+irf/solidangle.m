function angle=solidangle(VecA,VecB,VecC)
%IRF.SOLIDANGLE - calculates the solid angle with the sign taken into account.
%
%This function calculates the solid angle of three vectors making up a triangle in a unit sphere
%with the sign taken into account.
%
%   angle=irf.solidangle(VecA,VecB,VecC)
%
% Important: Time tags cannot be included in the input vectors 
% See Also C_4_POINCARE_INDEX
%Reference: R�e & Westergren (2004), Mathematics Handbook for Science and
%Engineering p.75 (Law of Cosine)

%--------written by E.Eriksson--------------------------------------------

%Check if time tags is in input vectors:
n=size(VecA,2); 
if n>3
error('Time tag not removed from the input vectors. Please do so and try again.')
end
if nargin==0
    help c_4_null_position;
    return;
elseif nargin<3
    error('Too few input values. Three vectors should be used.')
    elseif nargin>3
    error('Too many input values. Only three vectors should be used.')
end
%Calculate the smaller angles between the vectors around origin
a=acos(dot(VecC,VecB,2)); 
b=acos(dot(VecA,VecC,2)); 
c=acos(dot(VecB,VecA,2)); 
%Calculate the angles in the spherical triangle (Law of Cosines)
A=acos((cos(a)-cos(b).*cos(c))./(sin(b).*sin(c)));
B=acos((cos(b)-cos(a).*cos(c))./(sin(a).*sin(c)));
C=acos((cos(c)-cos(b).*cos(a))./(sin(b).*sin(a)));

%Calculates the Surface area on the unit sphere (solid angle)
angle=(A+B+C-pi); 
%Calculate the sign of the area
var=cross(VecC,VecB,2);
div=dot(var,VecA,2);
signarea=sign(div);

%Solid angle with sign taken into account
angle=signarea.*angle;
end