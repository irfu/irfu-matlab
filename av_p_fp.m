function [fp]=av_p_fp(n,m)
%function [fp]=av_p_fp(n,m)
% fp - plasma frequency in Hz
% n  - plasma density in cm^-3
% m  - mass in proton mass units, e.g.  m=1 (protons), m=2(Helium), m=0 (electrons)
% 

mp = 1.67e-27;
me = 9.11e-31;
e=1.6e-19;
eps=8.854e-12 ;
% fp=n e^2/eps/m/2/pi
fp=sqrt(e^2/eps)/2/pi;
fp=fp*sqrt(n*1e6);

if (m == 0)
   fp=fp/sqrt(me);
else
   fp=fp/sqrt(m*mp);
end

