function E2 = whistlerBtoE(B2,freq,thetak,Bmag,ne)
%WHISTLERBTOE compute electric field E power as a function of frequency for whistler
%waves using B power and cold plasma theory
%
% Input: 
%     B2 - power of whistler magnetic field in nT^2 Hz^{-1} (array or scalar)
%     freq - frequencies in Hz corresponding B2 (array or scalar)
%     thetak - wave-normal angle of whistler waves in degrees (scalar)
%     Bmag - magnetic field strength in nT (scalar)
%     ne - number density in cm^{-3} (scalar)
%
% Output:
%     E2 - electric field power (array or scalar)
%
% E.g. Epower = mms.whistlerBtoE(Bpower,freq,thetak,Bmag,ne)
% 
% Written by D. B. Graham

% Calculate plasma parameters
fpe = irf_plasma_calc(Bmag,ne,0,0,0,'Fpe');
fce = irf_plasma_calc(Bmag,ne,0,0,0,'Fce');

Units = irf_units;
c = Units.c;

%Check input
if length(B2) ~= length(freq)
    E2 = NaN;
    irf.log('critical','B2 and freq lengths do not agree!')
    return;
end

if size(B2) ~= size(freq)
    freq = freq';
end  

% Calculate cold plasma parameters
R = 1 - fpe^2./(freq.*(freq - fce));
L = 1 - fpe^2./(freq.*(freq + fce));
P = 1 - fpe^2./freq.^2;
D = 0.5*(R - L);
S = 0.5*(R + L);

n2 = R.*L*sind(thetak)^2 + P.*S*(1 + cosd(thetak)^2) - ...
    sqrt((R.*L - P.*S).^2*sind(thetak)^4 + 4*P.^2.*D.^2*cosd(thetak)^2);

n2 = n2./(2*(S*sind(thetak)^2 + P*cosd(thetak)^2));

Etemp1 = (P - n2*sind(thetak)^2).^2.*((D./(S - n2)).^2 + 1) + (n2*cosd(thetak)*sind(thetak)).^2;
Etemp2 = (D./(S - n2)).^2.*(P - n2*sind(thetak)^2).^2+P.^2*cosd(thetak)^2;

E2 = (c^2./n2).*Etemp1./Etemp2.*B2*1e-12;

end

