function res = irf_nv2p(n,v2)
%IRF_NV2P  Calculate plasma dynamic pressure
%
% res = irf_nv2p(n,v2)
%
% Calculate plasma dynamic pressure in nPa
% n in 1/cc
% v^2 in [km/s]^2
%
% p=nmv^2 ;-)

n = n(:);
v2 = v2(:);
% p=nmv^2 ;-)
res = 1.6726*1e-6*v2.*n;
end