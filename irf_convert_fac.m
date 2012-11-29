function [b_fac,e_fac]=irf_convert_fac(r,B,b,e)
%function [Rperpx,Rperpy,Rpar]=irf_convert_mfa(r,B,b,e) 
% 
% IRF_CONVERT_FAC transforms value A to a field-aligned coordinate
% (FAC) system
%   with R_parallel aligned with the background magnetic field
%   R_perp_y defined by R_parallel cross the position vector of the
%   spacecraft (nominally eastward at the equator)
%   R_perp_x defined by R_perp_y cross R_par
% r = position vector of spacecraft, columns (t x y z)
% B = background magnetic field, columns (t x y z)
% b = vector that is to be transformed to FAC, columns (t x y z)
% e = vector that is to be transformed to FAC, columns (t x y z)

% Note: all input parameters must be in the same coordinate system

  r=irf_resamp(r,B);
  b=irf_resamp(b,B);
  e=irf_resamp(e,B);
  

  %% Remove the last sample if the total number of samples is odd

  if size(r,1)/2 ~= floor(size(r,1)/2)
    r=r(1:end-1,:);
    B=B(1:end-1,:);
    b=b(1:end-1,:);
    e=e(1:end-1,:);
  end

  % set to zero NaNs
  ind_nan_r=isnan(r); r(ind_nan_r)=0;
  ind_nan_B=isnan(B); B(ind_nan_B)=0;
%  ind_nan_A=isnan(A); A(ind_nan_A)=0;
  
  
%% the direction of background magnetic field
bn=irf_norm(B);
t=r(:,1);
Rpar=bn;
Rperpy=irf_norm(irf_cross(Rpar, r));
Rperpx=irf_norm(irf_cross(Rperpy, B));

%A_mfa=A;
ndata=size(B,1);
b_fac=zeros(ndata,4);
b_fac(:,1)=b(:,1);
b_fac(:,4)=Rpar(:,2).*b(:,2)+Rpar(:,3).*b(:,3)+Rpar(:,4).*b(:,4);
b_fac(:,2)=Rperpx(:,2).*b(:,2)+Rperpx(:,3).*b(:,3)+Rperpx(:,4).*b(:,4);
b_fac(:,3)=Rperpy(:,2).*b(:,2)+Rperpy(:,3).*b(:,3)+Rperpy(:,4).*b(:,4);

e_fac=zeros(ndata,4);
e_fac(:,1)=e(:,1);
e_fac(:,4)=Rpar(:,2).*e(:,2)+Rpar(:,3).*e(:,3)+Rpar(:,4).*e(:,4);
e_fac(:,2)=Rperpx(:,2).*e(:,2)+Rperpx(:,3).*e(:,3)+Rperpx(:,4).*e(:,4);
e_fac(:,3)=Rperpy(:,2).*e(:,2)+Rperpy(:,3).*e(:,3)+Rperpy(:,4).*e(:,4);


end
