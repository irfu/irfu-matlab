% lp_sweep script to plot a standard IU curve for different antennas
%
% See also LP_PROBECURRENT


U_min=irf_ask('Umin(V) [%]>','U_min',-10); 
U_max=irf_ask('Umax(V) [%]>','U_max',20); 
Upot=U_min:0.1:U_max; Upot=Upot(:);
dU=0.01; % dU when estimating difference
Upot2=Upot+dU; ii=find(Upot<0);Upot2(ii)=Upot(ii)-dU;

R_sun=irf_ask('R_sun in AU [%]>','R_sun',0.22); 
m_amu1=irf_ask('mass of first species (in mp) [%]>','m_amu1',1); 
m_amu2=irf_ask('mass of second species (in mp) [%]>','m_amu2',16); 
m2=irf_ask('relative amount of second species  [%]>','m2',0); 
n_cc=irf_ask('Plasma density per cc [%]>','n_cc',100); 
n=n_cc*1e6; % n is plasma density per m3
Ti=irf_ask('Ion temperature in eV [%]>','Ti',20); 
Te=irf_ask('Electron temperature in eV [%]>','Te',20); 
V_SC=irf_ask('Spacecraft velocity in km/s [%]>','V_SC',0); 
UV_factor=irf_ask('UV factor [%]>','UV_factor',1); 
probe_type=irf_ask('Probe type. 1-spherical, 2- cylindrical 3-specify [%]>','probe_type',2); 
stazer_area=irf_ask('Probe total area in m2. [%]>','stazer_area',0.1885); % stazer surface area in m2
if probe_type==1,
    stazer_cross_section=stazer_area/pi;probe_text='spherical'; %sphere
    sunlit_total_surface_ratio=.25;
elseif probe_type ==2,
    stazer_cross_section=stazer_area/4;probe_text='cylindrical'; %cylinder
     sunlit_total_surface_ratio=1/pi;
elseif probe_type ==3,
    sunlit_total_surface_ratio=irf_ask('Specify ratio between sunlit and total area [%]:','sunlit_total_surface_ratio',.2);
    stazer_cross_section=stazer_area*sunlit_total_surface_ratio;probe_text='';probe_type=1; % for further calcualtions assume probe type sphere
end

J_probe=lp_probecurrent(probe_type,stazer_cross_section, ...
    stazer_area,Upot,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
J_probe_2=lp_probecurrent(probe_type,stazer_cross_section, ...
    stazer_area,Upot2,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
dUdI=(Upot2-Upot)./(J_probe_2-J_probe);

figure(71);
set(71,'Position',[10 10 600 1000]);
delete(findall(gcf,'Type','axes'))
h(1)=axes('position',[0.2 0.3 0.7 0.3]); % [x y dx dy]
h(2)=axes('position',[0.2 0.65 0.7 0.3]); % [x y dx dy]
h(3)=axes('position',[0.2 0.0 0.7 0.16]); % [x y dx dy]
info_txt='';
info_txt=[info_txt 'Ti=' num2str(Ti) 'eV, Te=' num2str(Te) 'eV, n=' num2str(n/1e6) 'cc, '];
info_txt=[info_txt 'R=' num2str(R_sun) 'AU.'];
info_txt=[info_txt '\newline ' probe_text ' probe area=' num2str(stazer_area) 'm2, area ratio sunlit/total=' num2str(sunlit_total_surface_ratio)];

hp=plot(h(1),Upot,J_probe*1e6);
grid(h(1),'on');
xlabel(h(1),'U [V]');
ylabel(h(1),'I [\mu A]');
set(h(1),'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
set(hp,'linewidth',2)

hp=plot(h(2),Upot,dUdI);
grid(h(2),'on');xlabel(h(2),'U [V]');
ylabel(h(2),'dU/dI [\Omega]');

Rmin = min(dUdI); % minimum resistance 
disp(['Minimum resistance                  R=' num2str(Rmin,3) ' Ohm']);
info_txt=[info_txt '\newline Minimum resistance R=' num2str(Rmin,3) ' Ohm'];
if min(J_probe)<0 && max(J_probe)>0,
   Ufloat=interp1(J_probe,Upot,0); % floating potential 
   info_txt=[info_txt '\newline Floating potential Ufl=' num2str(Ufloat,3) ' V'];
   disp(['Floating potential                Usc=' num2str(Ufloat,3) ' V']);
   Rfloat=interp1(Upot,dUdI,Ufloat);
   info_txt=[info_txt '\newline Resistance at floating potential R= ' num2str(Rfloat,3) ' Ohm'];
   disp(['Resistance at floating potential    R=' num2str(Rfloat,3) ' Ohm']);
end
set(h(2),'yscale','log')
set(h(2),'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
set(hp,'linewidth',2)

axis(h(3),'off');
text(0,1,info_txt,'fontsize',12);

