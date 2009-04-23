function res = irf_kev2k(T_eV,dir)
%IRF_KEV2K  convert Temperature from keV to K
%
%function res = irf_kev2k(T_eV,[dir])
%
% convert Temperature from keV to K
% dir=-1, konvert K to keV
% 
% See also: IRF_EV2MK
%
% $Id$

factor = 1.16*1e7;
if dir == -1
	factor = 1/factor;
end
if size(T_eV,2)==1
	res = factor*T_eV;
else
	res = T_eV;
	res(:,2) = factor*T_eV(:,2);
end