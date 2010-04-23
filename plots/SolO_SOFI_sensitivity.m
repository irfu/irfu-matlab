% sensitivity curve for the SolO antennas  
% based on Lennarts calculations for the Bepi


db_range=[-180 -40]; % dB range 
f_range=[1 1e6]; % frequencies in Hz
f_transition=[5e3 50e3]; % frequency range of transition from resistive to capacitive coupling

figure(1);clf;
h=subplot(1,1,1);
set(gca,'xscale','log','xlim',f_range);
set(gca,'ylim',db_range);
color_order=[0 0 1;1 0 0;0 0 1;0 0.6 0; 1 0 1];
set(gca,'linewidth',2,'MinorGridLineStyle','none','FontSize',14, ...
    'ColorOrder',color_order);
grid on;hold on;

preamp_noise=16e-9; % preamplifier noise 16nV/Hz1/2
antenna_eff_length=3; % efficient length of antenna from satellite center
preamp_noise_level=20*log10(preamp_noise/antenna_eff_length);
f_break=400; % transition frequency at which 1/f noise is starting
R_plasma_nobias=100e6; % 
R_plasma_bias=0.5e6; 
T_plasma=3*1e4; % 3eV
C_antenna=30e-12; % antenna capacitance in F
A_antenna=0.1885; % antenna area in m2
me= 9.1095e-031; % electron mass
n = 100e6; % plasma density
kB=1.3807e-023; % Boltzmann constant
qe=1.6022e-019; % charge 

% instrumental noise
SOFI_instr_noise_X=[f_range(1) f_break f_range(2)];
SOFI_instr_noise_Y=preamp_noise_level + [10*log10(f_break/f_range(1)) 0  0];

%plasma noise with bias
plasma_noise_bias=20*log10(sqrt(4*R_plasma_bias*T_plasma*kB));
SOFI_plasma_noise_bias_X=[f_range(1) f_transition(1)];
SOFI_plasma_noise_bias_Y=[plasma_noise_bias plasma_noise_bias];

%plasma noise without bias
plasma_noise_nobias=20*log10(sqrt(4*R_plasma_nobias*T_plasma*kB));
SOFI_plasma_noise_nobias_X=[f_range(1) f_transition(1)];
SOFI_plasma_noise_nobias_Y=[plasma_noise_nobias plasma_noise_nobias];

% shot noise 
nu=n/2*sqrt(8*kB*T_plasma/pi/me)*A_antenna;
f=10.^(log10(f_range(1)):.1:log10(f_range(2)));
SOFI_V_bias=10*log10(sqrt(2*nu*qe^2*R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length);
SOFI_V_nobias=10*log10(sqrt(2*nu*qe^2*R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length);

hp=plot(SOFI_instr_noise_X,SOFI_instr_noise_Y,'r',...
    f,SOFI_V_bias,'k', ...
    f,SOFI_V_nobias,'r', ...
    SOFI_plasma_noise_nobias_X,SOFI_plasma_noise_nobias_Y,'b', ...
    SOFI_plasma_noise_bias_X,SOFI_plasma_noise_bias_Y,'g');

set(hp,'linewidth',2);
ht=text(1e4,preamp_noise_level+5,'instrument');set(ht,'fontsize',14,'color','r');
ht=text(1e4,plasma_noise_bias+5,'biased');set(ht,'fontsize',14,'color','g');
ht=text(1e4,plasma_noise_nobias+5,'non-biased');set(ht,'fontsize',14,'color','b');

ylabel('Electric field intensity [dB V/m/Hz^{1/2}]');
xlabel('Frequency [Hz]');

patch([f_transition fliplr(f_transition)],db_range(1)+ [2 2 7 7], [0.8 0.8 0.8]);
patch([f_range(1) f_transition(1) f_transition(1) f_range(1)],db_range(1)+ [2 2 7 7], [0.9 0.9 0.9]);
ht=text(f_transition(1),db_range(1)+2,'resistive coupling    ');
set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','right');
patch([f_transition(2) f_range(2) f_range(2) f_transition(2)],db_range(1)+ [2 2 7 7], [0.9 0.9 0.9]);
ht=text(f_transition(2),db_range(1)+2,'   capacitive ');
set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left');
