function V=c_4_v(r1,r2,r3,r4,t)
%C_4_V  Calculate velocity V of discontinuity
%
%  v=C_4_V(r1,r2,r3,r4)  calculate velocity of discontinuity in
%  the same units per second as r
%
%  r1..r4 are vectors where first column is time when the satellite crosses
%         the discontinuity and the next three columns is satellite position
%
%  v=C_4_V(r1,r2,r3,r4,t)
%
%  r1..r4 are column vectors where first column is time and the next three
%         columns are satellite positions
%  t      is vector [t1 t2 t3 t4] where t1..t4 are times when satellite
%         cross the discontinuity
%
%  if more than 4 columns only the first 4 are used
%

%  Reference: ISSI book  Eq.12.9 Note that t_alpha should be(t_alpha-t_0)
%  m_l=Sk Rkl^-1
%  where Sk=(Sum dt? dr_?k)/N    N=4 number of satellites

if nargin<4;disp('Too few parameters. See usage:');help c_4_V;return;end
if nargin == 5
  r1=[t(1) interp1(r1(:,1),r1(:,[2 3 4]),t(1),'spline','extrap')];
  r2=[t(2) interp1(r2(:,1),r2(:,[2 3 4]),t(2),'spline','extrap')];
  r3=[t(3) interp1(r3(:,1),r3(:,[2 3 4]),t(3),'spline','extrap')];
  r4=[t(4) interp1(r4(:,1),r4(:,[2 3 4]),t(4),'spline','extrap')];
end

t1=r1(1,1);t2=r2(1,1);t3=r3(1,1);t4=r4(1,1);
t0=.25*(t1+t2+t3+t4);
dt1=t1-t0;dt2=t2-t0;dt3=t3-t0;dt4=t4-t0;

R1=r1(1,2:4);R2=r2(1,2:4);R3=r3(1,2:4);R4=r4(1,2:4);
R_Center=(R1+R2+R3+R4)/4;
dR1=R1-R_Center;dR2=R2-R_Center;dR3=R3-R_Center;dR4=R4-R_Center;

Rinv=c_4_r(dR1,dR2,dR3,dR4,-1); % inverse of volumetric tensor

S=(dt1*dR1+dt2*dR2+dt3*dR3+dt4*dR4)/4;

m=S*Rinv;

V=m/(m(1,1)^2+m(1,2)^2+m(1,3)^2);

if nargout==0
  disp(datestr(datenum(fromepoch(t(1)))))
  v=V;vn=irf_norm(v);dt=t-t(1);
  strv=['V=' num2str(irf_abs(v,1),3) '*[ ' num2str(vn(end-2:end),' %5.2f') '] km/s GSE'];
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  disp(strv);  disp(strdt);
end
