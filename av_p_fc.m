function [fc]=av_p_fc(B,m)
%function [fc]=av_p_fc(B,m)
% fc - cyclotron frequency in Hz
% B  - ambient field in nT
% mass in proton mass units, e.g. m=1 (protons), m=2(Helium), m=0 (electrons)
% 
mp = 1.67e-27;
me = 9.11e-31;
e=1.6e-19;

fc=e*B*1e-9/2/pi;

if (m == 0)
   fc=fc/me;
else
   fc=fc/m/mp;
end

