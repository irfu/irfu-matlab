function R=c_4_R(r1,r2,r3,r4,flag)
%c_4_R  Calculate volumetric tensor R
%  R=c_4_R(r1,r2,r3,r4,flag)  if flag=-1 returns the inverse of
%  volumetric tensor R, otherwise returns the volumetric tensor R.
%  flag can be omitted.
%  r1..r4 are row vectors with the satellite positions,
%  if more than 3 columns only last 3 are used
%
%  The volumetric tensor is calculated with respect to the center of r1..r4.
%  The units of volumetric tensor is the same as r

%  Reference: ISSI book
%  R=1/N Sum[sc=1->N] (r_sc -r_Center) (r_sc-r_Center)^T
%  r_Center = 1/N Sum[sc=1->N] r_sc
%
%

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'c_4_r')

if nargin<4;disp('See usage:');help c_4_R;return;end

if nargin>4,
   if flag ~= -1, flag=1;end  % calculate R if flag is not -1
elseif nargin<4,
   disp('ERROR! 4 s/c positions should be given. c_4_R()');
else
   flag=1; % calculate R if no flag given
end

R1=r1(1,end-2:end);R2=r2(1,end-2:end);R3=r3(1,end-2:end);R4=r4(1,end-2:end);

R_Center=(R1+R2+R3+R4)/4;
dR1=R1-R_Center;dR2=R2-R_Center;dR3=R3-R_Center;dR4=R4-R_Center;

R=(dR1'*dR1+dR2'*dR2+dR3'*dR3+dR4'*dR4)/4; % 1e6 to get into SI units km->m

if flag == -1,
R=inv(R);
end


