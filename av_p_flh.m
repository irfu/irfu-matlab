function [flh]=av_p_flh(B,n1,m1,n2,m2)
%function [flh]=av_p_flh(B,n1,m1)
%function [flh]=av_p_flh(B,n1,m1,n2,m2) % not implemented yet
% flh - lower hybrid frequency in Hz
% B   - ambient magnetic field in nT
% n1  - plasma density in cm^-3
% m1  - mass in proton mass units, e.g. m=1 (protons), m=2(Helium), m=0 (electrons)
% 
% See also AV_P_FP, AV_P_FC

mp = 1.67e-27;
me = 9.11e-31;
e=1.6e-19;
eps==8.854e-12 ;

% one component
flh=sqrt((av_p_fp(n1,m1)^2)*av_p_fc(B,0).^2./(av_p_fc(B,0).^2+av_p_fp(n1,0)^2)+av_p_fc(B,m1).^2);
% two component plasma
%av_p_fc(B,0)*2*pi*SQRT((D25+D33)/(H25+H26)+(D34*D33+D26*D25)/(D25+D33))/2/PI()

