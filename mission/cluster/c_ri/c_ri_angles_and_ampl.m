function [angles, ampl] = c_ri_angles_and_ampl(B1,B2,B3,B4)
%
% [angles, ampl] = c_ri_angles_and_ampl(B1,B2,B3,B4)
%
% Input: B1,B2,B3,B4 -[ time in epoch | Bx | By | Bz ]
%
% Output: angles -[ time in epoch |c1-c2|c1-c3|c1-c4|c2-c3|c2-c4|c3-c4]
%                 this are the angles between the cluster satellites
%         ampl   -[ampl c1 |ampl c2 |ampl c3 |ampl c4 ]
%                  the amplitude of the B-vectors
%
%Descrition of the function:
% calulates the angles between the cluster satellites and the amplitude
% of the cluster satellites.
%
%Using:
%
%Work method:
% formula: acos( M*N/(|M|*|N|) )*180/pi
%
%Error:
% NaN are put into the matrix where division by ZERO occurs.
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
warning off

if B1(1,1) ~= 0 && B2(1,1) ~= 0 && B3(1,1) ~= 0 && B4(1,1) ~= 0

  angles(:,1) = B1(:,1);

  [B1_n, a1] = norm_B(B1);
  [B2_n, a2] = norm_B(B2);
  [B3_n, a3] = norm_B(B3);
  [B4_n, a4] = norm_B(B4);

  ampl = [a1 a2 a3 a4];

  mDotn = M_dot_N(B1_n,B2_n);
  angles(:,2) = acos(mDotn);

  mDotn = M_dot_N(B1_n,B3_n);
  angles(:,3) = acos(mDotn);

  mDotn = M_dot_N(B1_n,B4_n);
  angles(:,4) = acos(mDotn);

  mDotn = M_dot_N(B2_n,B3_n);
  angles(:,5) = acos(mDotn);

  mDotn = M_dot_N(B2_n,B4_n);
  angles(:,6) = acos(mDotn);

  mDotn = M_dot_N(B3_n,B4_n);
  angles(:,7) = acos(mDotn);

  angles(:,2:7) = angles(:,2:7)*180/pi;

  angles(:,2:7) = real(angles(:,2:7));

else
  angles = 0;
  ampl = [0 0 0 0];
end

function [B_n, B_l] = norm_B(B)
B_l = sqrt((B(:,2).*B(:,2) + B(:,3).*B(:,3) + B(:,4).*B(:,4)));
B_m =[B_l B_l B_l];
B_n = B(:,2:4)./B_m;

%only B-value, no time
function B_dot = M_dot_N(M,N)
B_dot = (dot(M', N'))';

warning on
