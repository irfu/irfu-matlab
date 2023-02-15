function h=c_pl_sc_pos_mf(t,R1,R2,R3,R4,Vref,Bref) %#ok<INUSL>
%C_PL_SC_POS_MF   Plot spacecraft position in mean field coordinates
%
% that are defined usign spacecraft position and magnetic field.
% Z along magnetic field, Y axis is ZxR, X=YxZ where R is the
% position of spacecraft (s/c3)
%
% h=c_pl_sc_pos_mfinates(t,R1,R2,R3,R4,Vref,Bref);
% h=c_pl_sc_pos_mfinates(t);
% reads position, velocity and B field from mat files (use s/c 3 as reference)
%
% h    - handle to subplots
% t    - time in isdat epoch for which to make plot
% R1.. - s/c position (first column time, next columns x y z)
% V    - s/c velocity (first column time, next columns vx vy vz)
% B    - magnetic field (first column time, next columns Bx By Bz)
%

if nargin < 1, help c_pl_sc_pos_mf;return; end

figure;clf
if nargin == 1
  c_eval('load mR R?; r?=irf_resamp(R?,t,''spline'');clear R?;');
  c_load('V3');Vref=irf_resamp(V3,t);clear V3;
  c_load('B3');Bref=irf_resamp(B3,t);clear B3;
elseif nargin==7
  c_eval('r?=irf_resamp(R?,t,''spline'');clear R?;')
  Vref=irf_resamp(Vref,t);
  Bref=irf_resamp(Bref,t);
else
  help c_pl_sc_pos_mf;return;
end

%t=toepoch([2001 04 28 19 15 00]);
%load mBmod BT89Kp11; Bref=irf_resamp(BT89Kp11,t);clear BT89Kp11;

R=(r1+r2+r3+r4)/4;
c_eval('dr?=r?-R;dr?(1)=t;dr?=irf_abs(dr?);');
drref=max([dr1(5) dr2(5) dr3(5) dr4(5)]);
Vref_MF=irf_mean(Vref,r3,Bref);
Vscaling=.5*drref/irf_abs(Vref_MF,1);Vscale=Vscaling*Vref_MF(:,2:4);
c_eval('dr_MF?=irf_mean(dr?,r3,Bref);');

view1=[90 90]; %#ok<NASGU>
view2=[90 0]; %#ok<NASGU>
npl=2;ypl=2;xpl=1;
for ipl=1:npl
  h(ipl)=subplot(ypl,xpl,ipl);
  c_eval('x?=dr_MF?;');
  xl='outward [km] MF';yl='east [km] MF';zl='along B [km] MF';
  plot3(x1(2),x1(3),x1(4),'ks', ...
    x2(2),x2(3),x2(4),'rd', ...
    x3(2),x3(3),x3(4),'go', ...
    x4(2),x4(3),x4(4),'mv')
  xlabel(xl);ylabel(yl);zlabel(zl);
  c_eval('view(view?);',ipl);
  axis equal;grid on;
  line([0 Vscale(1)],[0 Vscale(2)],[0 Vscale(3)]);
end
axes(h(1));title(['v_{sc}=' num2str(irf_abs(Vref_MF,1),3) ' km/s, blue line r=v*t, t=' num2str(Vscaling,3) 's.'])



% subplot(ypl,xpl,1);
% cs_MF3=[-.49 .87];
% X=drref*[cs_MF3(1) -cs_MF3(1)];Y=drref*[cs_MF3(2) -cs_MF3(2)];
% hl=line(X,Y);set(hl,'Linewidth',2,'color','y');
% subplot(ypl,xpl,1);ax=1000;axis([-ax ax -ax ax -ax ax]);%axis equal;
% Y=[-500 -500-Vscaling*4];X=[-500 -500];
% hl=line(X,Y,[0 0]);set(hl,'color','b');
% subplot(ypl,xpl,2);ax=1000;axis([-ax ax -ax ax -ax ax]);%axis equal;
