clear hp hl;
NaturalConstants

Upot=-2:0.1:30; Upot=Upot(:);
dU=0.01; % dU when estimating difference
Upot2=Upot+dU; ii=find(Upot<0);Upot2(ii)=Upot(ii)-dU;

R_sun=1; % distance to sun in AU
m_amu1=1; % mass of first species in proton masses
m_amu2=16; % mass of second species in proton masses
m2=0; % relative number of second species

there_is_more_curves=1;j=0;
while there_is_more_curves
  j=j+1;
  n_cc=irf_ask('Plasma density per cc [%]>','n_cc',100);
  n=n_cc*1e6; % n is plasma density per m3
  Ti=irf_ask('Ion temperature in eV [%]>','Ti',20);
  Te=irf_ask('Electron temperature in eV [%]>','Te',20);
  V_SC=0; % probe velocity with respect to media
  UV_factor=1; % default is 1
  probe_type=irf_ask('Probe type. 1-spherical, 2- cylindrical [%]>','probe_type',2);
  probe_radius=irf_ask('Probe radius in cm. [%]>','probe_radius',4);
  probe_area=4*pi*(probe_radius/100)^2;
  if probe_type==1
    probe_cross_section=probe_area/pi; %cylinder
  elseif probe_type ==2
    probe_cross_section=probe_area/4; %sphere
  end


  J_probe=lp_probecurrent(probe_type,probe_cross_section, ...
    probe_area,Upot,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
  J_probe_2=lp_probecurrent(probe_type,probe_cross_section, ...
    probe_area,Upot2,R_sun,UV_factor,m_amu1,m_amu2,m2,n,Ti,Te,V_SC);
  dUdI=(Upot2-Upot)./(J_probe_2-J_probe);
  J_probe=J_probe*1e6;% to get in microA
  J_probe_2=J_probe_2*1e6;% to get in microA

  result.Ti(j)=Ti;
  result.Te(j)=Te;
  result.n(j)=n;
  result.Upot{j}=Upot;
  result.dUdI{j}=dUdI;
  result.Jprobe{j}=J_probe;

  there_is_more_curves=irf_ask('There is more UI curves to plot? (1/0) [%]>','there_is_more_curves',0);

end

number_of_curves=j;

figure;
%title_txt=['Ti=' num2str(Ti) 'eV ' 'Te=' num2str(Te) 'eV '  ' probe radius=' num2str(probe_radius) 'cm'];
%h(1)=irf_subplot(2,1,-1);
ccc=get(gca,'colororder');
ccc=[0 0 0;ccc]; % add balck color first
for j=1:number_of_curves
  hp(j)=plot(result.Upot{j},result.Jprobe{j},'color',ccc(j,:));
  grid on;hold on;
  legend_txt{j}=['n=' num2str(result.n(j)/1e6) 'cc ' 'Ti=' num2str(result.Ti(j)) 'eV ' 'Te=' num2str(result.Te(j)) 'eV'];
end

set(gca,'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
set(gca,'ylim',[-0.3999 0.099]);
set(gca,'xlim',[-1.99 29.99]);
set(hp,'linewidth',2)

switch number_of_curves
  case 1
    legend(legend_txt{1},'location','southeast');
  case 2
    legend(legend_txt{1},legend_txt{2},'location','southeast');
  case 3
    legend(legend_txt{1},legend_txt{2},legend_txt{3},'location','southeast');
  case 4
    legend(legend_txt{1},legend_txt{2},legend_txt{3},legend_txt{4},'location','southeast');
end
set(hl,'fontsize',12);

ylabel('I [\mu A/m2]');
xlabel('U [V]');
ht=irf_pl_info([' probe radius ' num2str(probe_radius) 'cm'],gca,[0,.92]);
set(ht,'FontSize', 12);
ht=irf_pl_info('Xscale_IU.m',gca,[0,1]);
set(ht,'interpreter','none','FontSize', 8);

add_ref_points=irf_ask('Add reference points? y/n>','qq','n');
while strcmp(add_ref_points,'y')
  iref=irf_ask('Reference points at current level I[micro A/m2]=','ilevel',0);
  for j=1:number_of_curves
    uref=interp1(result.Jprobe{j},result.Upot{j},iref);
    plot(uref,iref,'o','color',ccc(j,:),'MarkerFaceColor',ccc(j,:),'MarkerSize',8);
  end
  add_ref_points=irf_ask('Additional reference points? y/n [%]>','add_ref_points','n');
end

%ht=text(10,-0.07,'unbiased probes');set(ht,'fontsize',12);
%ht=text(4,-0.2,'biased probes'); set(ht,'fontsize',12)

% in case resistance also needed
% h(2)=irf_subplot(2,1,-2);
% for j=1:number_of_curves,
%     hp(j)=plot(result.Upot{j},result.dUdI{j},'color',ccc(j,:));
%     grid on;hold on;
%     disp(['Minimum resistance R=' num2str(min(result.dUdI{j}),3) ' Ohm']);
% end
%
% xlabel('U [V]');ylabel('R [Ohm]');
% set(h(2),'yscale','log')
% set(h(2),'linewidth',2,'MinorGridLineStyle','none','FontSize',14)
% set(hp,'linewidth',2)
%
%print -depsc2 IU_so_T20_n100_A0.35.eps
