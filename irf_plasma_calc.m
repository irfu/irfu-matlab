function [Wpe,Wce,Wuh,Wpp,Wcp,WpO,WcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,noshow)
%IRF_PLASMA_CALC   Calculate basic plasma quantities
%
% irf_plasma_calc(B,n,no,Te,Ti)
%
% [Wpe,Wce,Wuh,Wpp,Wcp,WpO,WcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,flag);
%
%	B - magnetic field [nT]
%	n - density [cm^-3]
%	no - content of O+ [%]
%	Te - electron temperature [eV]
%	Ti - ion temperature [eV]
% flag - 1: do not display the output
%
% $Id$

% Copyright 1997-2004 Yuri Khotyaintsev

if nargin < 6
	noshow = 0;
end

Mp_Me = 1836.15;

Wpe = 8.973*sqrt(n);
Wce = 2.8e-2*B;
Wpp = 8.973*sqrt(n*(1-no/100)/Mp_Me);  
WpO = 8.973*sqrt(n*(no/100)/Mp_Me/16); 
Va = 21.8*B/sqrt(n*(1+15*no/100));
Vte = 4.19*1e2*sqrt(Te);
Vtp = 9.79*sqrt(Ti);
Vts = 9.79*sqrt(Te);
VtO = Vtp/4;
Le = 5.3e3/sqrt(n);
Li = 3e8/(2*pi*Wpp*1e3);

Wpe = Wpe*1e3;
Wce = Wce*1e3;
Wuh = sqrt(Wce^2+Wpe^2);
Wpp = Wpp*1e3;
Wcp = Wce/Mp_Me;
WpO = WpO*1e3;
WcO = Wce/Mp_Me/16;

Rop = Vtp/Wcp*1e3; % in meters
RoO = VtO/WcO*1e3;
Ros = Vts/Wcp*1e3;

if noshow > 0
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies

freq = [Wpe Wce Wuh Wpp Wcp WpO WcO];
freqs = {'W_pe'; 'W_ce'; 'W_uh'; 'W_pp'; 'W_cp'; 'W_pO'; 'W_cO'};
disp(sprintf('\nPlasma frequencies\n'))
for ii = 1:length(freq)
	val = freq(ii);
	if val > 1e6
		units = '[MHz]';
		koef = 1e-6;
	elseif val < 1e6 & val > 1e3
		units = '[KHz]';
		koef = 1e-3;
    elseif val <1e3 & val >.1
		units = '[Hz]';
		koef = 1;
    else
        units = '[mHz]';
		koef = 1e3;
	end
	disp(sprintf('%s = %5.2f %s',freqs{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velosities

v = [Va Vte Vtp Vts VtO ];
vs = {'V_a'; 'V_Te'; 'V_Tp'; 'C_s'; 'V_TO'};

disp(sprintf('\nPlasma velosities\n'))
for ii = 1:length(v)
	val = v(ii);
	if val > 1e3
		units = '[km/s 1e3]';
		koef = 1e-3;
	elseif val < 1e3 & val > 1
		units = '[km/s]';
		koef = 1;
	else
		units = '[m/s]';
		koef = 1e3;
	end
	disp(sprintf('%s = %5.2f %s',vs{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scales

l = [Le Li Rop RoO Ros];
ls = {'L_e'; 'L_i'; 'Ro_p'; 'Ro_O'; 'Ro_s'};
disp(sprintf('\nPlasma scales\n'))
for ii = 1:length(l)
	val = l(ii);
	if val > 1e3
		units = '[km]';
		koef = 1e-3;
	elseif val < 1e3 & val > 1
		units = '[m]';
		koef = 1;
	else
		units = '[cm]';
		koef = 1e2;
	end
	disp(sprintf('%s = %5.2f %s',ls{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% etc

disp(sprintf('\nPlasma dimensionless parameters\n'))

beta  = Vtp^2/Va^2;

if beta < 10/Mp_Me
	disp(sprintf('beta*Mp_Me = %1.5f',Vtp^2/Va^2*Mp_Me))
else
	disp(sprintf('beta  = %1.5f',Vtp^2/Va^2))
end
