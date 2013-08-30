%
% A matlab implementation of the PR2000.FOR
% returns the precession matrix 
% Ex:
% p=pr2000(jd2k)
%
%
% 20061127, Verified and confirmed with fortran routines
%
function p=pr2000(mjd2k)

T = mjd2k - 0.5D0;

GZ = T*(0.3061153D-6 + T*(0.10976D-14 + T*0.179D-20));
ZA = GZ + T*T*(0.2881D-14 + T*0.358D-22);
TH = T*(0.2660417D-6 - T*(0.1550D-14 + T*0.41549D-20));


CGZ=cos(GZ);
SGZ=sin(GZ);
CZA=cos(ZA);
SZA=sin(ZA);
CTH=cos(TH);
STH=sin(TH);
    
p=[CGZ*CZA*CTH - SGZ*SZA, -SGZ*CZA*CTH - CGZ*SZA,-CZA*STH;
   CGZ*SZA*CTH + SGZ*CZA, -SGZ*SZA*CTH + CGZ*CZA,-SZA*STH;
   CGZ*STH,-SGZ*STH,CTH];
