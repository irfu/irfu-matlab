% sensitivity curve for the SolO antennas  
% based on Lennarts calculations for the Bepi

db_range=[-180 -80]; % dB range 
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

bias_noise_level=-180+20*log10(30/5); % noise level  
SOFI_bias_X=[f_range(1) 1e3 f_range(2)];
SOFI_bias_Y=bias_noise_level + [log10(SOFI_bias_X(2)/SOFI_bias_X(1))*10 0  0];

resistance_ratio=10; % how many times unbiased resistance is larger than biased
nobias_noise_level=-180+20*log10(30/5)+20*log10(resistance_ratio); % noise level
SOFI_nobias_X=[f_range(1) 1e3 f_range(2)];
SOFI_nobias_Y=nobias_noise_level + [log10(SOFI_nobias_X(2)/SOFI_nobias_X(1))*10 0  0];

hp=plot(SOFI_nobias_X,SOFI_nobias_Y,'r',SOFI_bias_X,SOFI_bias_Y,'b');
set(hp,'linewidth',2);
ht=text(1e4,bias_noise_level+5,'biased');set(ht,'fontsize',14,'color','b');
ht=text(1e4,nobias_noise_level+5,'non-biased');set(ht,'fontsize',14,'color','r');

ylabel('Electric field intensity [dB V/m/Hz^{1/2}]');
xlabel('Frequency [Hz]');

patch([f_transition fliplr(f_transition)],db_range(1)+ [2 2 7 7], [0.8 0.8 0.8]);
patch([f_range(1) f_transition(1) f_transition(1) f_range(1)],db_range(1)+ [2 2 7 7], [0.9 0.9 0.9]);
ht=text(f_transition(1),db_range(1)+2,'resistive coupling    ');
set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','right');
patch([f_transition(2) f_range(2) f_range(2) f_transition(2)],db_range(1)+ [2 2 7 7], [0.9 0.9 0.9]);
ht=text(f_transition(2),db_range(1)+2,'   capacitive ');
set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left');



%title('Solar Orbiter RPW noise level');

% for jj=1:length(v);
%     ttt=['v= ' num2str(v_kms(jj)) ' km/s'];
%     ht=text(l_inv_range(1)*11,f_range(2)*0.95*0.5^jj,ttt);
%     set(ht,'FontSize',14,'VerticalAlignment','top','Color',color_order(jj,:));
% end

% RPW instrument DC;LF;HF
% patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [f_range(1) f_range(1) 300 300],[1 0.5 0.5]) 
% ht=text(l_inv_range(2)*0.12, 10,'DC');set(ht,'fontsize',14);
% patch(l_inv_range(2)*[0.33 0.99 0.99 0.33], [100 100 25000 25000],[1 0.5 0.5]) 
% ht=text(l_inv_range(2)*0.4, 1000,'LF');set(ht,'fontsize',14);
% patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [2.5e3 2.5e3 1e6 1e6],[1 0.5 0.5]) 
% ht=text(l_inv_range(2)*0.12, 1e5,'HF');set(ht,'fontsize',14);

