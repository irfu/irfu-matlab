function res = irf_ev2mk(T_eV,dir)
%IRF_EV2MK  convert Temperature from eV to MK
%
%function res = irf_ev2mk(T_eV,[dir])
%
% convert Temperature from keV to K
% dir=-1, konvert MK to eV
%
% See also: IRF_KEV2K
%
% $Id$


disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp(' Will be removed')
disp(' Use irf_convert')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')


factor = 1.1604*1e-2;
if ~( nargin==1 || dir ~= -1)
	factor = 1/factor;
end
if size(T_eV,2)==1
	res = factor*T_eV;
else
	res = T_eV;
	res(:,2) = factor*T_eV(:,2);
end