function h=c_plot_sc_position_mf_coord(t,R1,R2,R3,R4,Vref,Bref);
% C_PLOT_SC_POSITION_MF_COORD - plot spacecraft position in mean field coordinates
% that are defined usign spacecraft position and magnetic field.
% Z along magnetic field, Y axis is ZxR, X=YxZ where R is the position of spacecraft (s/c3)
%
% h=c_plot_sc_position_mf_coordinates(t,R1,R2,R3,R4,Vref,Bref);
% h=c_plot_sc_position_mf_coordinates(t);   % reads position, velocity and B field from mat files (use s/c 3 as reference)
%
% h    - handle to subplots
% t    - time in isdat epoch for which to make plot
% R1.. - s/c position (first column time, next columns x y z)
% V    - s/c velocity (first column time, next columns vx vy vz)
% B    - magnetic field (first column time, next columns Bx By Bz)
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'c_pl_sc_pos_mf')

if nargin < 1, help c_plot_sc_position_mf_coord;return; end

figure;clf
if nargin == 1,
  for ic=1:4,eval(av_ssub('load mR R?; r?=av_interp(R?,t,''spline'');clear R?;',ic)),end
  load mV V3;Vref=av_interp(V3,t);clear V3;
  load mB B3;Bref=av_interp(B3,t);clear B3;
elseif nargin==7,
  for ic=1:4,eval(av_ssub('r?=av_interp(R?,t,''spline'');clear R?;',ic)),end
  Vref=av_interp(Vref,t);
  Bref=av_interp(Bref,t);
else
  help c_plot_sc_position_mf_coord;return;
end

%t=toepoch([2001 04 28 19 15 00]);
%load mBmod BT89Kp11; Bref=av_interp(BT89Kp11,t);clear BT89Kp11;

R=(r1+r2+r3+r4)/4;
for ic=1:4,eval(av_ssub('dr?=r?-R;dr?(1)=t;dr?=av_abs(dr?);',ic)),end
drref=max([dr1(5) dr2(5) dr3(5) dr4(5)]);
%Vref_MF=av_mean(Vref,r3,Bref,'GSE');
Vref_MF=av_mean(Vref,r3,Bref);
Vscaling=.5*drref/av_abs(Vref_MF,1);Vscale=Vscaling*Vref_MF(:,2:4);
% for ic=1:4,eval(av_ssub('dr_MF?=av_mean(dr?,r3,Bref,''GSE'');',ic)),end
for ic=1:4,eval(av_ssub('dr_MF?=av_mean(dr?,r3,Bref);',ic)),end

view1=[90 90];
view2=[90 0];
npl=2;ypl=2;xpl=1;
for ipl=1:npl,
  h(ipl)=subplot(ypl,xpl,ipl);
  for ic=1:4,eval(av_ssub('x?=dr_MF?;',ic)),end;xl='outward [km] MF';yl='east [km] MF';zl='along B [km] MF';
  plot3(x1(2),x1(3),x1(4),'ks', ...
       x2(2),x2(3),x2(4),'rd', ...
        x3(2),x3(3),x3(4),'go', ...
        x4(2),x4(3),x4(4),'mv')
  xlabel(xl);ylabel(yl);zlabel(zl);
  eval(av_ssub('view(view?);',ipl));
  axis equal;grid on;
  line([0 Vscale(1)],[0 Vscale(2)],[0 Vscale(3)]);
end
axes(h(1));title(['v_{sc}=' num2str(av_abs(Vref_MF,1),3) ' km/s, blue line r=v*t, t=' num2str(Vscaling,3) 's.'])



% subplot(ypl,xpl,1);
% cs_MF3=[-.49 .87];
% X=drref*[cs_MF3(1) -cs_MF3(1)];Y=drref*[cs_MF3(2) -cs_MF3(2)];
% hl=line(X,Y);set(hl,'Linewidth',2,'color','y');
% subplot(ypl,xpl,1);ax=1000;axis([-ax ax -ax ax -ax ax]);%axis equal;
% Y=[-500 -500-Vscaling*4];X=[-500 -500];
% hl=line(X,Y,[0 0]);set(hl,'color','b');
% subplot(ypl,xpl,2);ax=1000;axis([-ax ax -ax ax -ax ax]);%axis equal;
