function [Blnm1,h]=p1052_ar(tint,Vmp,L,Nsign,flagEdotB);
% function [Blnm]=pl_cluster_ejxb(tint,Vmp,L,Nsign,flagEdotB);
%
% Magnetic field in LMN
% electric field in comparisoin to jxB
% tangential electric field
% satellite potential
% current from single s/c method
%
% Vmp - magnetopause velocity, N=norm(Vmp)
% L - L vector
% Nsign - 1 if N in the direction of Vmp, -1 if N in the direction 
%         opposite to Vmp (e.g., N shall always point from MSph to Msh)
% flagEdotB - 1 if use E.B=o assumption for the 3rd component
%
% Usage: 
%   [blnm] = pl_cluster_ejxb(toepoch([2002 03 30 13 11 42])+[0 8],...
%          31.4* [ -0.94 -0.21 -0.25],[0 0 1],-1,0);
%
% $Id$

persistent Vmp_str L_str tint_str sc_list icb  Te 

c_load('B?');
c_load('diB?');
c_load('diE?p1234');
c_load Damp?
c_load Ddsi?
c_eval('diE?p1234 = caa_corof_dsi(diE?p1234,Ddsi?,Damp?);')
%load mE dvE1 dvE2 dvE3 dvE4;
c_load('P?');
%load mP P1 P2 P3 P4 NVps1 NVps2 NVps3 NVps4;

if nargin <1,
  help p1052;
  tint_str=irf_ask('Time interval [%]>','tint_str','toepoch([2002 03 30 13 11 42])+[0 8]');
  tint=eval(tint_str);
end
if nargin <2,
  Vmp_str=irf_ask('Magnetopause velocty Vmp (N=Nsign*norm(Vmp)) [%]>','Vmp_str','1 0 0');
  Vmp=eval(['[' Vmp_str ']']);
end
if nargin <3,
  L_str=irf_ask('Input L vector, (M=NxL,L=MxN) [%]>','L_str','0 0 1');
  L=eval(['[' L_str ']']);
end
if nargin <4,
  Nsign=irf_ask('Input Nsign [%]>','Nsign',1);
end

N=Nsign*Vmp;
M=cross(L,N);L=cross(N,M);
N=N./norm(N);M=M./norm(M);L=L./norm(L);

flag_EdotB=1; % assume E.B=0 to calculate spin axis E component
v=-Vmp; % velocity of sc with respect to discontinuity
dt =c_v([tint(1) Vmp]); dt = dt - dt(2);

tint_data=tint+[min(dt) max(dt)];
c_eval('B?=irf_tlim(irf_abs(B?),tint_data);');
c_eval('diB?=irf_tlim(diB?,tint_data);');
c_eval('dvE?=irf_tlim(diE?p1234,tint_data);');

if flag_EdotB==1,
  c_eval('[dvE?,d?]=irf_edb(dvE?,diB?,10);');
end

c_eval('dEvxb?=irf_tappl(irf_cross(irf_resamp(diB?,dvE?),c_gse2dsc([dvE?(1,1) Vmp],?,2)),''*1e-3*(-1)'');');
c_eval('dvE?=irf_add(1,dvE?,1,dEvxb?);');

c_eval('dNsp?=c_gse2dsc([B?(1,1) N],?,2);dNsp?(1,4)=0;dNsp?=irf_norm(dNsp?);dMsp?=irf_cross([dNsp?(1,1) 0 0 -1],dNsp?);');  % the direction in spin plane closest to the magnetopause normal
c_eval('dN?=c_gse2dsc([B?(1,1) N],?,2);dN?=irf_norm(dN?);');
c_eval('dM?=c_gse2dsc([B?(1,1) M],?,2);dM?=irf_norm(dM?);');
c_eval('dEn?=irf_dot(dvE?,dNsp?);dEm?=irf_dot(dvE?,dMsp?);');  % the direction in spin plane closest to the magnetopause normal
if flag_EdotB==1,
  c_eval('dEn?=irf_dot(dvE?,dN?);dEm?=irf_dot(dvE?,dM?);');
end

c_eval('[j?]=irf_jz(v,B?);jz?=irf_dot(j?,irf_norm(B?));jtot?=irf_abs(j?,1);jperp?=jz?;jperp?(:,2)=sqrt(jtot?.^2-jz?(:,2).^2);n?=c_efw_scp2ne_20040104(P?);');

c_eval('ejb?=irf_vec_x_scal(irf_tappl(irf_cross(j?,B?),''*1e-9*1e3''),irf_tappl(n?,''*1.6e-19*1e6''),-1);'); % j=[A],B=[nT],n=[cc],E=[mV/m]
c_eval('dejb?=c_gse2dsc(ejb?,?,2);');
c_eval('dejbn?=irf_dot(dejb?,dNsp?);');

c_eval('Blnm?=irf_newxyz(B?,L,N,M);');

%Te=irf_ask('Electron energy Te?[%]>','Te',1);
Te1=3000; %PS
Te2=50;   %FTE
Ti1=6000; %PS
Ti2=300;  %FTE

% Remove bad points due to spin
for ax=3:4
	if ax==3, tx = toepoch([2004 1 4 12 47 7.2]);
	else, tx = toepoch([2004 1 4 12 47 7.6]);
	end
	c_eval('ii=find(n?(:,1)>tx & n?(:,1)<tx+.6);n?(ii,:)=[];',ax)
end
  
c_eval('dtn?=n?;dtn?(1,2)=0;dtn?(2:end,1)=0.5*(dtn?(2:end,1)+dtn?(1:end-1,1));');
c_eval('dtn?(2:end,2)=(n?(2:end,2)-n?(1:end-1,2))./(n?(2:end,1)-n?(1:end-1,1));');
c_eval('vn?=irf_resamp([0 v],dtn?);');
% Approved formula for -grad(P)/ne = T/(ne)*dtn*Nsign/abs(Vmp)
c_eval('egradn?=[dtn?(:,1) Nsign*Te1./irf_abs(vn?,1)./n?(:,2).*dtn?(:,2)];');
c_eval('igradn?=[dtn?(:,1) -Nsign*Ti1./irf_abs(vn?,1)./n?(:,2).*dtn?(:,2)];');
c_eval('egradn1?=[dtn?(:,1) Nsign*Te2./irf_abs(vn?,1)./n?(:,2).*dtn?(:,2)];');
c_eval('igradn1?=[dtn?(:,1) -Nsign*Ti2./irf_abs(vn?,1)./n?(:,2).*dtn?(:,2)];');
%c_eval('egradn3d?=irf_vec_x_scal(irf_tappl(irf_norm(vn?),''*(-1)''),egradn?);');

sc_list=irf_ask('For which s/c plot single s/c plots ?[%]>','sc_list',[]);

icb=irf_ask('For which s/c to include internal burst ?[%]>','icb',[]);

for j=1:length(icb),
  sc_id=icb(j);
  c_load('P4kHz?p4',sc_id);
  c_eval('nburst?=c_efw_scp2ne_20040104(P4kHz?p4);',sc_id);
  %c_eval('nburst?(1:2,:)=[];nburst?(:,2)=lowpass(nburst?(:,2),30,9000);',sc_id);
  %c_eval('dtn?=nburst?;dtn?(1,2)=0;dtn?(2:end,1)=0.5*(dtn?(2:end,1)+dtn?(1:end-1,1));',sc_id);
  %c_eval('dtn?(2:end,2)=(nburst?(2:end,2)-nburst?(1:end-1,2))./(nburst?(2:end,1)-nburst?(1:end-1,1));',sc_id);
  %c_eval('vn?=irf_resamp([0 v],dtn?);',sc_id);
  % Approved formula for -grad(P)/ne
  %c_eval('egradn_burst?=[dtn?(:,1) Nsign*Te./irf_abs(vn?,1)./nburst?(:,2).*dtn?(:,2)];',sc_id);
  %c_eval('egradn3d_burst?=irf_vec_x_scal(irf_tappl(irf_norm(vn?),''*(-1)''),egradn_burst?);',sc_id);
  c_eval('t0=nburst?(1,1);te=nburst?(end,1);ntmp=nburst?;',sc_id);
  c_eval('ii = find(n?(:,1)<t0); ntmp = [ n?(ii,:); ntmp];',sc_id);
  c_eval('n? = ntmp; clear ntmp',sc_id);
end



%%%%%%%%%%%%%%  Figure 1 %%%%%%%%%%%%%%
%%%%%%%%%%%%%% All sc %%%%%%%%%%%%%%%%%

figure;h=irf_plot({Blnm1,Blnm1,Blnm2,Blnm3,Blnm3,Blnm3,Blnm3,Blnm3,Blnm3});

axes(h(1));
c_pl_tx('Blnm?',2,dt.*0);ylabel('B_L [nT]');
ht=irf_pl_info([mfilename '  ' datestr(now) ...
  ' Vmp=' num2str(irf_abs(Vmp,1)) ' [' num2str(Vmp(1:3),' %6.2f') ']km/s,' ...
  ' dt=[' num2str(dt,' %6.2f') ']s,' ...
  ' L=[' num2str(L,' %6.2f') '], N=[' num2str(N,' %6.2f') '].'],gca,[0,1 ]);
leg_coord=[.1,.6];font_size=13;cluster_labels={'C1 ','C2 ','C3 ','C4 '};cluster_colors='krgb';
for ic=1:4,
  ht=irf_pl_info(cluster_labels{ic},gca,leg_coord);set(ht,'color',cluster_colors(ic),'Fontsize',font_size,'FontWeight','demi');
  ext=get(ht,'extent'); leg_coord=leg_coord+[ext(3)*1.2 0];
end

axes(h(2));
c_pl_tx('Blnm?',2,dt);ylabel('B_L [nT]');

axes(h(3));
c_pl_tx('Blnm?',4,dt);ylabel('B_M [nT]');

axes(h(4));
c_pl_tx('Blnm?',3,dt);ylabel('B_N [nT]');

axes(h(5));
c_pl_tx('irf_abs(Blnm?)',5,dt);ylabel('B_{tot} [nT]');

axes(h(6));
cluster_colors='krgb';
c_eval('irf_plot(n?,''dt'',dt(?),''LineStyle'',cluster_colors(?));hold on;',[1 3 4 2])
%c_pl_tx('n?',2,dt);
ylabel('N_{Vps} [cc]');
set(h(6),'YScale','log')
set(h(6),'YLim',[.3 19])

axes(h(7));
c_pl_tx(dEm2(1,:),dEm2,dEm3,dEm4,2,dt);
ylabel('E_{M} [mV/m]');
%sht=irf_pl_info(['f=20-180Hz'],gca,[0.05,.8]); set(ht,'interpreter','none','FontSize', 10);

axes(h(8));
c_pl_tx('dEn?',2,dt);
ylabel('E_{N} [mV/m]');
%ylabel('E,jxB/ne \newline [mV/m]');
%hold on;
%c_pl_tx('dejbn?',2,dt,{'.','.','.','.'});

axes(h(9))
c_eval('fs?=c_efw_fsample(dEn?);')
dd1=0;dd2=0;dd3=.5;dd4=-.5;
c_eval('eint?=[dEn?(:,1) 1e-3*irf_abs(Vmp,1)/fs?*cumsum(dEn?(:,2))+dd?];')
c_pl_tx('eint?',2,dt);
ylabel('U [kV]')
set(h(9),'YLim',[0 5.9],'YTick',[1 2 3 4 5])

lw=1;k=-3:0;
for j=1:9,
  hh=get(h(j),'children');set(hh(end+[k k]),'linewidth',lw);
end

axis(h,'tight');
irf_zoom(tint,'x',h);
add_timeaxis(h);

irf_pl_mark([toepoch([2004 01 04 12 47 10]) + [-.15 .4]],h);
irf_pl_mark([toepoch([2004 01 04 12 47 8.4]) + [-.05 .05]],h);

%ht=irf_pl_number_subplots(h);

%%%%%%%%%%%%%%  Figure single s/c %%%%%%%%%%%%%%
%%%%%%%%%%%%%% Separate s/c %%%%%%%%%%%%%%%%%

for ic=sc_list;
  figure;
  h=irf_plot({Blnm1,Blnm1,Blnm2,Blnm3});

  axes(h(1));
  %if find(ic==icb), % plot also burst data
  %  c_eval('irf_plot({dEn?,dejbn?,egradn?,egradn_burst?},''comp''); ylabel(''E,jxB/ne,\nabla p/ne\newline [mV/m] sc?'');',ic);
  %else
    c_eval('irf_plot({dEn?,dejbn?,egradn?},''comp''); ylabel(''E,jxB/ne,\nabla p/ne \newline[mV/m] sc?'');',ic);
  %end
  %c_eval('irf_plot(Blnm?(:,[1 5])); ylabel(''B [nT] sc?''); ',ic);
  ht=irf_pl_info([mfilename '  ' datestr(now) ...
    ' Vmp=[' num2str(Vmp(1:3),' %6.2f') ']km/s,' ...
    ' dt=[' num2str(dt,' %6.2f') ']s,' ...
    ' L=[' num2str(L,' %6.2f') '], N=[' num2str(N,' %6.2f') '].'],gca,[0,1 ]);
  leg_coord=[.2,.6];font_size=13;cluster_labels={'C1 ','C2 ','C3 ','C4 '};cluster_colors='krgb';

  axes(h(2));
  c_eval('irf_plot(n?(:,[1 2]));',ic);
  if find(ic==icb), % plot also burst data
    c_eval('hold on; irf_plot(nburst?);',ic);
  end
  c_eval('ylabel(''N_{Vps} [cc] sc?'');',ic);
  set(h(2),'YScale','log')

  axes(h(3));
  c_eval('irf_plot(Blnm?(:,[1 5])); ylabel(''B [nT] sc?''); ',ic);

  axes(h(4));
  c_eval('irf_plot(irf_tappl([jz? jperp?(:,2)],''*1e6''));ylabel(''j_{||},jperp [\mu A/m^2] sc?'');',ic);

  axis(h,'tight');
  irf_zoom(tint,'x',h);
  add_timeaxis(h);

  %ht=irf_pl_number_subplots(h)

end
%%%%%%%%%%%%%%  Figure 3 %%%%%%%%%%%%%%
%%%%%%%%%%%%%% Separate s/c %%%%%%%%%%%%%%%%%



figure;
h=irf_plot({Blnm1,Blnm1,Blnm2,Blnm3,Blnm3,Blnm3});

axes(h(1));
%c_pl_tx('n?',2,dt);
cluster_colors='krgb';
c_eval('irf_plot(n?,''dt'',dt(?),''LineStyle'',cluster_colors(?));hold on;',[1 3 4 2])
hold off
ylabel('N_{Vps} [cc]');
set(h(1),'YScale','log')
set(h(1),'YLim',[.3 19])
ht=irf_pl_info([mfilename '  ' datestr(now) ...
  ' Vmp=[' num2str(Vmp(1:3),' %6.2f') ']km/s,' ...
  ' dt=[' num2str(dt,' %6.2f') ']s,' ...
  ' L=[' num2str(L,' %6.2f') '], N=[' num2str(N,' %6.2f') '].'],gca,[0,1 ]);
leg_coord=[.1,.6];font_size=13;cluster_labels={'C1 ','C2 ','C3 ','C4 '};
for ic=1:4,
  ht=irf_pl_info(cluster_labels{ic},gca,leg_coord);set(ht,'color',cluster_colors(ic),'Fontsize',font_size,'FontWeight','demi');
  ext=get(ht,'extent'); leg_coord=leg_coord+[ext(3)*1.2 0];
end

axes(h(2))
c_eval('jz? = irf_tappl(jz?,''*1e6'');')
c_pl_tx('jz?',2,dt);
ylabel('j_{||}[\mu A/m^2]');
leg_coord=[.1,.6];font_size=13;cluster_labels={'C1 ','C2 ','C3 ','C4 '};cluster_colors='krgb';
for ic=1:4,
  ht=irf_pl_info(cluster_labels{ic},gca,leg_coord);set(ht,'color',cluster_colors(ic),'Fontsize',font_size,'FontWeight','demi');
  ext=get(ht,'extent'); leg_coord=leg_coord+[ext(3)*1.2 0];
end
set(h(2),'YTick',[-.4 -.2 0 .2 .4])

for ax=1:4,
  axes(h(ax+2));
  tmax = toepoch([2004 1 4 12 47 9.6]);
  c_eval('ii=find(egradn?(:,1)>tmax+dt(?));egradn?(ii,:)=[];')
  c_eval('ii=find(igradn?(:,1)>tmax+dt(?));igradn?(ii,:)=[];')
  tmax = toepoch([2004 1 4 12 47 9.8]);
  c_eval('ii=find(egradn1?(:,1)<tmax+dt(?));egradn1?(ii,:)=[];')
  c_eval('ii=find(igradn1?(:,1)<tmax+dt(?));igradn1?(ii,:)=[];')
  c_eval(['irf_plot(dEn?,''dt'',dt(?));hold on;'...
  'irf_plot(dejbn?,''dt'',dt(?),''LineStyle'',''g'');hold on;'...
  'irf_plot(egradn?,''dt'',dt(?),''LineStyle'',''r'');'...
  'irf_plot(egradn1?,''dt'',dt(?),''LineStyle'',''r'');'...
  'irf_plot([egradn?(end,:); egradn1?(1,:)],''dt'',dt(?),''LineStyle'',''r--'');'...
  'irf_plot(igradn?,''dt'',dt(?),''LineStyle'',''k'');'...
  'irf_plot(igradn1?,''dt'',dt(?),''LineStyle'',''k'');'...
  'irf_plot([igradn?(end,:); igradn1?(1,:)],''dt'',dt(?),''LineStyle'',''k--'');'...
  'ylabel(''[mV/m] C?'');'],ax);
  hold off
  set(gca,'YLim',[-39 64])
  
  leg_coord=[.1,.6];font_size=13;cluster_labels={'E ', 'uxB ','jxB/ne','-\nabla p/ne'};cluster_colors='bkgr';
  for ic=1:4,
    ht=irf_pl_info(cluster_labels{ic},gca,leg_coord);set(ht,'color',cluster_colors(ic),'Fontsize',font_size,'FontWeight','demi');
    ext=get(ht,'extent'); leg_coord=leg_coord+[ext(3)*1.2 0];
  end
end

%axis(h,'tight');
%irf_zoom(tint,'x',h);
irf_zoom(toepoch([2004 01 04 12 47 07]) + [.5 4],'x',h)
add_timeaxis(h);

lw=1;k=-5:0;
for j=1:2
  hh=get(h(j),'children');
  set(hh(end+[k k]),'linewidth',lw);

end
lw=1;k=-7:0;
for j=3:6
  hh=get(h(j),'children');
  set(hh(end+[k k]),'linewidth',lw);
  set(gca,'YLim',[-39 64])
end

irf_pl_mark([toepoch([2004 01 04 12 47 10]) + [-.12 .4]],h);
irf_pl_mark([toepoch([2004 01 04 12 47 8.4]) + [-.05 .05]],h);
%ht=irf_pl_number_subplots(h)
%irf_zoom(toepoch([2004 01 04 12 47 10]) + [-2.5 1.2],'x',h);
