NaturalConstants

Upot=-10:0.1:20; Upot=Upot(:);
dU=0.01; % dU when estimating difference
Upot2=Upot+dU; ii=find(Upot<0);Upot2(ii)=Upot(ii)-dU;

R_sun=.22, % distance from sun in AU
m_amu1=1; % mass of first species in proton masses
m_amu2=16; % mass of second species in proton masses
m2=0; % relative number of second species
n=100e6; % plasma density per m3
Ti=20; % ion temperature in eV
Te=20; % ion temperature in eV
V_SC=0; % probe velocity with respect to media
UV_factor=1; % default is 1
probe_type=2; % cylindrical 
stazer_area=0.35; % stazer surface area in m2
stazer_cross_section=stazer_area/pi; % 
J_probe=lp_probecurrent(probe_type,stazer_cross_section, ...
    stazer_area,Upot,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
J_probe_2=lp_probecurrent(probe_type,stazer_cross_section, ...
    stazer_area,Upot2,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
dUdI=(Upot2-Upot)./(J_probe_2-J_probe);
J_probe=J_probe*1e6;% to get in microA
J_probe_2=J_probe_2*1e6;% to get in microA

figure;
title_txt=['Ti=' num2str(Ti) 'eV, Te=' num2str(Te) 'eV, n=' num2str(n/1e6) 'cc, R=' num2str(R_sun) 'AU. cylindrical probe area=' num2str(stazer_area) 'm2'];
h(1)=subplot(2,1,1);
hp=plot(Upot,J_probe);grid on;xlabel('U [V]');ylabel('I [\mu A/m2]');
title(title_txt)
set(gca,'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
set(hp,'linewidth',2)

h(2)=subplot(2,1,2);
hp=plot(Upot,dUdI);grid on;xlabel('U [V]');ylabel('dUdI [Ohm]');
set(h(2),'yscale','log')
set(h(2),'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
set(hp,'linewidth',2)
%print -depsc2 IU_so_T20_n100_A0.35.eps
