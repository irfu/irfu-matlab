function [out]=irf_mean(inp,r,b,z)
% function [out]=irf_mean(inp,r,b)
% function [out]=irf_mean(inp,r,b,z)
%
% put inp into mean field coordinates defined by
% position vector r and magnetic field b
% if earth magnetic dipole axis z is given then uses another algorithm (good for auroral passages)
%
% inp = Nx4 vector, first column time
% r = Nx4 vector, first column time
% b = Nx4 vector, first column time
% z = Nx4 vector, first column time or 'GSE' (it means use earth dipole in that coordinate)
% out = Nx4 vector, first column time
% if r,z,b of different lengths than inp then do interpolation of r,b,z
% Z axis is along b, Y axis is Zxr, X=YxZ
% if z is given then
% Z axis is along b, Y axis is zxb.sign(b.r), X=YxZ

out=inp;flag_dipole=0;
if nargin > 3
 flag_dipole=1;
 if size(z,1) ~= size(inp,1)
  zz=irf_resamp(z,inp);
 else
  zz=z;
 end
end
if (size(r,1) ~= size(inp,1))
    rr=irf_resamp(r,inp);
else
    rr=r;
end
if (size(b,1) ~= size(inp,1))
    bb=irf_resamp(b,inp);
else
    bb=b;
end
zv=irf_norm(bb);
if flag_dipole == 0
 yv=irf_norm(irf_cross(zv,rr));
else
 ss=irf_dot(b,r);ind= ss(:,2) > 0;
 ss(:,2)=-1;ss(ind,2)=1;
 yv=irf_norm(irf_vec_x_scal(irf_cross(zz,bb),ss));
end
xv=irf_cross(yv,zv);
%keyboard
% in case rotation axis is used as reference uncomment next line
% rot_axis=rr;rot_axis(:,[2 3])=0;yv=irf_norm(irf_cross(irf_cross(bb,rot_axis),bb));xv=irf_cross(yv,zv);

out(:,2)=irf_dot(xv,inp,1);
out(:,3)=irf_dot(yv,inp,1);
out(:,4)=irf_dot(zv,inp,1);
