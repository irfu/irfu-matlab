function [div,grad]=c_4_gradvec(r1,r2,r3,r4)
%c_4_gradVec  Calculate spatial gradient tensor of vector (div and curl as special cases)
%  [gradvec]=c_4_gradvec(r1,r2,r3,r4)  calculate spatial gradient tensor of vector
%  the gradient is in the same units as r
%
%  [div,curl]=c_4_gradvec(r1,r2,r3,r4)  calculate divergence and curl of vector
%
%  r1..r4 are row vectors where
%  column 1     is time
%  column 2-4   are satellite position
%  column 5-7   are vector values
%
%  Reference: ISSI book  Eq.12.18
%  k_l=S_l R_kl^-1
%  where S_l=(Sum a!=b dx_ab dR_ab)/N^2   N=4 number of satellites
%

if nargin<4;    disp('Too few parameters. See usage:');help c_4_gradVec;     return;end

R1=r1(1,2:4);R2=r2(1,2:4);R3=r3(1,2:4);R4=r4(1,2:4);
R_Center=(R1+R2+R3+R4)/4;
dR1=R1-R_Center;dR2=R2-R_Center;dR3=R3-R_Center;dR4=R4-R_Center;
x1=r1(1,5:7);x2=r2(1,5:7);x3=r3(1,5:7);x4=r4(1,5:7);

S=0;
for a=1:4, for b=1:a,     S=S+eval(av_ssub('(x?-x!)''*(dR?-dR!)',a,b));      end,end
S=S/16;

Rinv=c_4_R(dR1,dR2,dR3,dR4,-1);                                     % inverse of volumetric tensor

gradVec=S*Rinv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0&nargout==1,
   div=gradVec;
elseif nargout == 2,
       div=trace(gradVec);
       curl=[gradVec(2,3)-gradVec(3,2) -gradVec(1,3)+gradVec(3,1) gradVec(1,2)-gradVec(2,1)];
end

