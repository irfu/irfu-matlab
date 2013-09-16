function res = irf_e_vxb(v,b,flag)
%IRF_E_VXB   Compute VxB and ExB/B^2
%
% [E] = IRF_E_VXB(V,B)
%        Calculate electric field E = -VxB
%
% [VExB] = IRF_E_VXB(E,B,-1)
%     Calculate convection velocity VExB = ExB/B^2
%
%     Units: v[km/s], B[nT], E[mV/m]
%     V = [T VX VY VZ] (T..VZ column vectors)
%     B = [T BX BY BZ]
%     V and B can be at different sampling
%     Resulting E is at B sampling, VExB is at E sampling
%
% $Id$

if (nargin ==3) && (flag == -1)
	% Estimate VExB = ExB/B^2
	e = v;
	if size(b,1) ~= size(e,1), b = irf_resamp(b,e); end
	res = irf_vec_x_scal( irf_cross(e,b), [b(:,1) irf_abs(b,1)], -2 );
	res = irf_tappl(res,'*1e3');
else
	% Estimate E =(v x B)
	if size(v,1) == 1
		if size(v,2) == 3
			v = [b(1,1) v];
		elseif size(v,2) < 3
			error('v has to few components');
		end
		v = irf_resamp(v,b);
	else
		b = irf_resamp(b,v);
	end
	res = irf_tappl( irf_cross(v,b), '*1e3*1e-9*1e3*(-1)' );
end

