function [p,z,f,fim,P1] = m2xyz(m)
% M2XYZ orders the output from WHAMP into matrices
%
% 	[p,z,f,fim,P1] = WHAMP.M2XYZ(m)
%
% m  - matrix where the first and second column are wave vector components p and z
%       normally the third column would be freqyency and
%       the fourth would be imaginary part of frequency and later whatever
% p  - perpendicular wave vector component
% z  - parallel wave vector component
% f  - matrix with real part of frequencies
% fim - matrix with imaginary part of frequencies
% P1 - [optional] matrix of parameter in the 5th column
%
% You can adopt routine to even more parameters if needed 
%

if nargin<1, help whamp.m2xyz;return;end

n = length(m);
width = length(m(1,:));
number_output = min(width,nargout);
ix = 1; p(1)=m(1,1);
iy = 1; z(1)=m(1,2);
for i = 1:n,
 if m(i,1) ~= p(ix)
  if m(i,1) ~= p(1)
   ix = ix +1;
   p(ix) = m(i,1);
  else
   ix = 1;
  end
 end
 if m(i,2) ~= z(iy)
  if m(i,2) ~= z(1)
   iy = iy +1;
   z(iy) = m(i,2);
  else
   iy = 1;
  end
 end
 f(iy,ix) = m(i,3);
 if (number_output > 3)
  	fim(iy,ix) = m(i,4);
  	if number_output > 4
  		P1(iy,ix)=m(i,4);
  	end
 end
end
