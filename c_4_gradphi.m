function gradPhi=c_4_gradPhi(r1,r2,r3,r4)
%c_4_gradPhi  Calculate spatial gradient of scalar
%  gradPhi=c_4_gradPhi(r1,r2,r3,r4)  calculate spatial gradient of scalar
%  the gradient is in the same units as r
%  r1..r4 are row vectors where
%  column 1     is time
%  column 2-4   are satellite position
%  column 5     is value of scalar
%
%  Reference: ISSI book  Eq.12.16
%  k_l=S_l R_kl^-1
%  where S_l=(Sum a!=b dx_ab dRab)/N^2   N=4 number of satellites
%

if nargin<4;disp('Too few parameters. See usage:');help c_4_gradPhi;return;end

R1=r1(1,2:4);R2=r2(1,2:4);R3=r3(1,2:4);R4=r4(1,2:4);
R_Center=(R1+R2+R3+R4)/4;
dR1=R1-R_Center;dR2=R2-R_Center;dR3=R3-R_Center;dR4=R4-R_Center;
x1=r1(1,5);x2=r2(1,5);x3=r3(1,5);x4=r4(1,5);

S=0;
for a=1:4, for b=1:a,     S=S+eval(av_ssub('(x?-x!)*(dR?-dR!)',a,b));      end,end
S=S/16;

Rinv=c_4_R(dR1,dR2,dR3,dR4,-1);                                     % inverse of volumetric tensor

gradPhi=S*Rinv;

