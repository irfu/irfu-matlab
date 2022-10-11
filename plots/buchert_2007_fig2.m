function [h]=buchert_2007_fig2(tint,Vmp,L,Nsign,flagEdotB)
% function [h]=buchert_2007_fig2(tint,Vmp,L,Nsign,flagEdotB);
%
% based on p1052.m
%
% Magnetic field in LMN
% electric field in comparisoin to jxB
% tangential electric field
% satellite potential
% current from single s/c method
%
% Vmp - magnetopause velocity, N=norm(Vmp)
% L - L vector
% Nsign - 1 if N in the direction of Vmp, -1 if N in the direction opposite to Vmp (e.g., N shall always point from MSph to Msh)
% flagEdotB - 1 if use E.B=o assumption for the 3rd component
%
% Usage: buchert_2007_fig2
% complicated untested alternative: [blnm]=buchert_2007_fig2(toepoch([2002 03 30 13 11 42])+[0 8],31.4* [ -0.94 -0.21 -0.25],[0 0 1],-1,0);

persistent Vmp_str L_str tint_str ib ic

c_load('B?');
c_load('diB?');
c_load('diE?p1234');
%load mE dvE1 dvE2 dvE3 dvE4;
c_load('P?');
c_load('NVps?');
%load mP P1 P2 P3 P4 NVps1 NVps2 NVps3 NVps4;

if nargin <1
  help buchert_2007_fig2;
  tint_str=irf_ask('Time interval [%]>','tint_str','toepoch([2002 03 04 09 48 34])+[0 6]');
  tint=eval(tint_str);
end
if nargin <2
  Vmp_str=irf_ask('Magnetopause velocty Vmp (N=Nsign*norm(Vmp)) [%]>','Vmp_str','-90.3 -38.7 13.9');
  Vmp=eval(['[' Vmp_str ']']);
end
if nargin <3
  L_str=irf_ask('Input L vector, (M=NxL,L=MxN) [%]>','L_str','0.13 0.06 0.99');
  L=eval(['[' L_str ']']);
end
if nargin <4
  Nsign=irf_ask('Input Nsign [%]>','Nsign',1);
end

N=Nsign*Vmp;
M=cross(L,N);L=cross(N,M);
N=N./norm(N);M=M./norm(M);L=L./norm(L);

flag_EdotB=1; % assume E.B=0 to calculate spin axis E component
v=-Vmp; % velocity of sc with respect to discontinuity
dt =c_v([tint(1) Vmp]);

tint_data=tint+[min(dt) max(dt)]; %#ok<NASGU>
c_eval('B?=irf_tlim(irf_abs(B?),tint_data);');
c_eval('diB?=irf_tlim(diB?,tint_data);');
c_eval('dvE?=irf_tlim(diE?p1234,tint_data);');
c_eval('NVps?=irf_tlim(NVps?,tint_data);');

if flag_EdotB==1
  c_eval('[dvE?,d?]=irf_edb(dvE?,diB?,10);');
end

% callibrate dvE? and subtract MP motion
% offs=0.5*[1 1 1 1];coef=1.5* [1 1 1 1];
%c_eval('dvE?(:,2)=dvE?(:,2)+offs(?);dvE?(:,2:end)=dvE?(:,2:end)*coef(?);');

c_eval('dEvxb?=irf_tappl(irf_cross(irf_resamp(diB?,dvE?),c_gse2dsc([dvE?(1,1) Vmp],?,2)),''*1e-3*(-1)'');');
c_eval('dvE?=irf_add(1,dvE?,1,dEvxb?);');

c_eval('dNsp?=c_gse2dsc([B?(1,1) N],?,2);dNsp?(1,4)=0;dNsp?=irf_norm(dNsp?);dMsp?=irf_cross([dNsp?(1,1) 0 0 -1],dNsp?);');  % the direction in spin plane closest to the magnetopause normal
c_eval('dN?=c_gse2dsc([B?(1,1) N],?,2);dN?=irf_norm(dN?);');
c_eval('dM?=c_gse2dsc([B?(1,1) M],?,2);dM?=irf_norm(dM?);');
c_eval('dEn?=irf_dot(dvE?,dNsp?);dEm?=irf_dot(dvE?,dMsp?);');  % the direction in spin plane closest to the magnetopause normal
if flag_EdotB==1
  c_eval('dEn?=irf_dot(dvE?,dN?);dEm?=irf_dot(dvE?,dM?);');
end

c_eval('[j?]=irf_jz(v,B?);jz?=irf_dot(j?,irf_norm(B?));n?=c_efw_scp2ne(P?);');
c_eval('ejb?=irf_vec_x_scal(irf_tappl(irf_cross(j?,B?),''*1e-9*1e3''),irf_tappl(n?,''*1.6e-19*1e6''),-1);'); % j=[A],B=[nT],n=[cc],E=[mV/m]
c_eval('dejb?=c_gse2dsc(ejb?,?,2);');
c_eval('dejbn?=irf_dot(dejb?,dNsp?);');


c_eval('Blnm?=irf_newxyz(B?,L,N,M);');

%%%%%%%%%%%%%%  Figure 2 %%%%%%%%%%%%%%
%%%%%%%%%%%%%% Separate s/c %%%%%%%%%%%%%%%%%
ic=irf_ask('which s/c to plot?[%]>','ic',2);
Te=irf_ask('Electron energy Te?[%]>','Te',100);

c_eval('dtn=NVps?;dtn(1,2)=0;dtn(2:end,1)=0.5*(dtn(2:end,1)+dtn(1:end-1,1));',ic);
c_eval('dtn(2:end,2)=(NVps?(2:end,2)-NVps?(1:end-1,2))./(NVps?(2:end,1)-NVps?(1:end-1,1));',ic);
c_eval('vn=irf_resamp([0 v],dtn);',ic);
c_eval('egradn=[dtn(:,1) Te./irf_abs(vn,1)./NVps?(:,2).*dtn(:,2)];',ic);
c_eval('egradn3d=irf_vec_x_scal(irf_tappl(irf_norm(vn),''*(-1)''),egradn);',ic);

figure;
h=irf_plot({Blnm1,Blnm1,Blnm2,Blnm3});

axes(h(1));cla
ht=irf_pl_info([mfilename '  ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss")) ...
  ' Vmp=' num2str(irf_abs(Vmp,1)) ' [' num2str(Vmp(1:3),' %6.2f') ']km/s,' ...
  ' dt=[' num2str(dt,' %6.2f') ']s,' ...
  ' L=[' num2str(L,' %6.2f') '], N=[' num2str(N,' %6.2f') '].'],gca,[0,1 ]);
leg_coord=[.2,.6];font_size=13;cluster_labels={'C1 ','C2 ','C3 ','C4 '};cluster_colors='krgb';

axes(h(2));
c_eval('irf_plot(Blnm?(:,[1 5])); ylabel(''B [nT] sc?''); ',ic);

axes(h(3));
if strcmp(ib,'y')
  c_eval('irf_plot({dEn?,dejbn?,egradn,egradn_burst},''comp''); ylabel(''E,jxB/ne,\nabla p/ne [mV/m] sc?'');',ic);
else
  c_eval('irf_plot({dEn?,dejbn?,egradn},''comp''); ylabel(''E,jxB/ne,\nabla p/ne [mV/m] sc?'');',ic);
end
irf_zoom([-12 7],'y')

axes(h(4));
c_eval('irf_plot([irf_tappl(jz?,''*1e6'') irf_abs(j?,1).*1e6]);ylabel(''|j|,jz [\mu A/m2] sc?'');',ic);

axis(h,'tight');
irf_zoom(tint,'x',h);
irf_timeaxis(h);

numb={'a','b','c','d','e','f','g','h','i','j','k','l','mf'};
for ip=1:4
  axes(h(ip));
  ht=irf_pl_info(numb{ip},gca,[0.01,.8]);
  set(ht,'fontsize',12);
end
