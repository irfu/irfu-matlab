function [b_fac,e_fac]=irf_convert_fac(r,B0,b,e)
%IRF_CONVERT_FAC transforms to a field-aligned coordinate (FAC) system
%
%  [b_fac,e_fac] = irf_convert_mfa(r,B0,b,[e])
%
%  Transforms to a field-aligned coordinate (FAC) system defined as:
%  R_parallel_z aligned with the background magnetic field
%  R_perp_y defined by R_parallel cross the position vector of the
%  spacecraft (nominally eastward at the equator)
%  R_perp_x defined by R_perp_y cross R_par
%
%  Input:
%    r  = position vector of spacecraft, columns (t x y z)
%    B0 = background magnetic field, columns (t x y z)
%    b  = vector that is to be transformed to FAC, columns (t x y z)
%    e  = vector that is to be transformed to FAC, columns (t x y z)

% Note: all input parameters must be in the same coordinate system

if nargin<4, e=[];
elseif ~isempty(e) && (size(b,1) ~= size(e,1))
    error('E and B must be of the same size')
end
if size(b,1) ~= size(B0,1), B0 = irf_resamp(B0,b); end
if size(b,1) ~= size(B0,1), r = irf_resamp(r,b); end  
  
  
%% the direction of background magnetic field
bn=irf_norm(B0);
Rpar=bn;
Rperpy=irf_norm(irf_cross(Rpar, r));
Rperpx=irf_norm(irf_cross(Rperpy, B0));

%A_mfa=A;
ndata=size(B0,1);
b_fac=zeros(ndata,4);
b_fac(:,1)=b(:,1);
b_fac(:,4)=Rpar(:,2).*b(:,2)+Rpar(:,3).*b(:,3)+Rpar(:,4).*b(:,4);
b_fac(:,2)=Rperpx(:,2).*b(:,2)+Rperpx(:,3).*b(:,3)+Rperpx(:,4).*b(:,4);
b_fac(:,3)=Rperpy(:,2).*b(:,2)+Rperpy(:,3).*b(:,3)+Rperpy(:,4).*b(:,4);

if ~isempty(e)
    e_fac=zeros(ndata,4);
    e_fac(:,1)=e(:,1);
    e_fac(:,4)=Rpar(:,2).*e(:,2)+Rpar(:,3).*e(:,3)+Rpar(:,4).*e(:,4);
    e_fac(:,2)=Rperpx(:,2).*e(:,2)+Rperpx(:,3).*e(:,3)+Rperpx(:,4).*e(:,4);
    e_fac(:,3)=Rperpy(:,2).*e(:,2)+Rperpy(:,3).*e(:,3)+Rperpy(:,4).*e(:,4);
end

end
