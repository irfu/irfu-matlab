function [varargout]=mms4_gradvec(r1,r2,r3,r4,u1,u2,u3,u4)
%C_4_GRADVEC calculate spatial gradient tensor of vector
%
%  [gradvec]=mms4_gradvec(R1,R2,R3,R4,x1,x2,x3,x4)  calculate spatial gradient 
%  tensor of vector. The gradient is in the same units as r.
%
%  [div,curl]=mms4_gradvec(R1,R2,R3,R4,x1,x2,x3,x4)  calculate divergence and curl of vector
%
%  r1..r4 SC position in TSeries format
%  u1..u4 vector in TSeries format
%
%  Reference: ISSI book  Eq.12.18
%  k_l=S_l R_kl^-1
%  where S_l=(Sum a!=b dx_ab dR_ab)/N^2   N=4 number of satellites
%
% $Id$

if nargin<4;    disp('Too few parameters. See usage:');help c_4_gradVec;     return;end
    


n=size(u1.data,1);
gradVec_data=zeros(n,3,3);
for i=1:n % to be vectorized at some point
    
R1=r1.data(i,1:3);R2=r2.data(i,1:3);R3=r3.data(i,1:3);R4=r4.data(i,1:3);
R_Center=(R1+R2+R3+R4)/4;
dR1=R1-R_Center;dR2=R2-R_Center;dR3=R3-R_Center;dR4=R4-R_Center;
x1=u1.data(i,1:3);x2=u2.data(i,1:3);x3=u3.data(i,1:3);x4=u4.data(i,1:3);

S=0;
for a=1:4, for b=1:a,     S=S+eval(irf_ssub('(x?-x!)''*(dR?-dR!)',a,b));      end,end
S=S/16;

Rinv=c_4_r(dR1,dR2,dR3,dR4,-1);                                     % inverse of volumetric tensor

gradVec_data(i,:,:)=S*Rinv;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
   gradVec = TSeries(u1.time,gradVec_data,'TensorOrder',2,...
           'repres',{'x','y','z'},'repres',{'x','y','z'});


if nargout==0||nargout==1,
varargout(1)={gradVec};
elseif nargout == 2,
       div=trace(gradVec);
       curl=[gradVec(2,3)-gradVec(3,2) -gradVec(1,3)+gradVec(3,1) gradVec(1,2)-gradVec(2,1)];
       varargout(1)={div};
       varargout(2)={curl};
end
