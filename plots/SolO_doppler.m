% plot figure showing typical temporal and spatial scales and doppler
% shifts for solar orbiter plasma conditions

NaturalConstants

f_range=[0.1 1e6]; % frequencies in Hz
l_range=[1 1e6]; % wave length in m 
l_inv_range=fliplr(1./l_range); % inverse wave length
l_deb_range=[2 10]; l_inv_deb_range=fliplr(1./l_deb_range);
l_erho_range=[100 500]; l_inv_erho_range=fliplr(1./l_erho_range);
l_e_range=[200 2000]; l_inv_e_range=fliplr(1./l_e_range);
l_rho_range=[4e3 40e3]; l_inv_rho_range=fliplr(1./l_rho_range);
l_h_range=[7e3 70e3]; l_inv_h_range=fliplr(1./l_h_range);
f_cp_range=[.07 7]; % H+ gyrofrequency
f_ce_range=[140 14e3]; % H+ gyrofrequency
f_pe_range=[30e3 300e3]; % H+ gyrofrequency
f_lh_range=[3 300]; % lower hybrid freq
v_kms=[300 1000 1e4 3e5]; % velocity lines in km/s
v=v_kms*1e3; % to get in m (SI units)

figure(1);clf;
h=subplot(1,1,1);
set(gca,'xscale','log','xlim',l_inv_range);
set(gca,'yscale','log','ylim',f_range);
color_order=[0 0 0;1 0 0;0 0 1;0 0.6 0; 1 0 1];
set(gca,'linewidth',2,'MinorGridLineStyle','none','FontSize',14, ...
    'ColorOrder',color_order);
grid on;hold on;
f_v_doppler_lines=repmat(l_inv_range,length(v),1);
f_v_doppler_lines(:,1)=f_v_doppler_lines(:,1).*v';
f_v_doppler_lines(:,2)=f_v_doppler_lines(:,2).*v';
hp=plot(l_inv_range,f_v_doppler_lines);
set(hp,'linewidth',2);
ylabel('f [Hz]');
xlabel('1/\lambda [m^{-1}]');
for jj=1:length(v);
    ttt=['v= ' num2str(v_kms(jj)) ' km/s'];
    ht=text(l_inv_range(1)*11,f_range(2)*0.95*0.5^jj,ttt);
    set(ht,'FontSize',14,'VerticalAlignment','top','Color',color_order(jj,:));
end

% RPW instrument DC;LF;HF
patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [f_range(1) f_range(1) 300 300],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.12, 10,'DC');set(ht,'fontsize',14);
patch(l_inv_range(2)*[0.33 0.99 0.99 0.33], [100 100 25000 25000],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.4, 1000,'LF');set(ht,'fontsize',14);
patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [2.5e3 2.5e3 1e6 1e6],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.12, 1e5,'HF');set(ht,'fontsize',14);

% proton inertial length
patch([l_inv_h_range fliplr(l_inv_h_range)], f_range(1)*[3.2 3.2 8 8],[0.5 1 0.5]) 
ht=text(l_inv_h_range(1)*1.12,f_range(1)*5,'c/\omega_{pi}');set(ht,'fontsize',14);
% proton gyror
patch([l_inv_rho_range fliplr(l_inv_rho_range)], f_range(1)*[1.25 1.25 3 3],[0.5 1 0.5]) 
ht=text(l_inv_rho_range(1)*1.12,f_range(1)*2,'\rho_H+');set(ht,'fontsize',14);
% electron inertial
patch([l_inv_e_range fliplr(l_inv_e_range)], f_range(1)*[3.2 3.2 8 8],[0.5 1 0.5]) 
ht=text(l_inv_e_range(1)*1.12,f_range(1)*5,'c/\omega_{pe}');set(ht,'fontsize',14);
% electron gyror
patch([l_inv_erho_range fliplr(l_inv_erho_range)], f_range(1)*[1.25 1.25 3 3],[0.5 1 0.5]) 
ht=text(l_inv_erho_range(1)*1.12,f_range(1)*2,'\rho_e');set(ht,'fontsize',14);
% debye length
patch([l_inv_deb_range fliplr(l_inv_deb_range)], f_range(1)*[1.25 1.25 3 3],[0.5 1 0.5]) 
ht=text(l_inv_deb_range(1)*1.12,f_range(1)*2,'\lambda_D');set(ht,'fontsize',14);

% H+ gyrofreq
patch( l_inv_range(1)*[1.2 1.2 2.6 2.6], [f_cp_range fliplr(f_cp_range)],[0.8 0.8 1]) 
ht=text(l_inv_range(1)*1.23,10^(sum(log10(f_cp_range))/2),'f_{cH+}');set(ht,'fontsize',14);
% lower hybrid freq
patch( l_inv_range(1)*[2.8 2.8 6 6], [f_lh_range fliplr(f_lh_range)],[0.8 0.8 1]) 
ht=text(l_inv_range(1)*3,10^(sum(log10(f_lh_range))/2),'f_{LH}');set(ht,'fontsize',14);
% e- gyrofreq
patch( l_inv_range(1)*[1.2 1.2 2.6 2.6], [f_ce_range fliplr(f_ce_range)],[0.8 0.8 1]) 
ht=text(l_inv_range(1)*1.23,10^(sum(log10(f_ce_range))/2),'f_{ce}');set(ht,'fontsize',14);
% e- plasma freq.
patch( l_inv_range(1)*[2.8 2.8 6 6], [f_pe_range fliplr(f_pe_range)],[0.8 0.8 1]) 
ht=text(l_inv_range(1)*3,10^(sum(log10(f_pe_range))/2),'f_{pe}');set(ht,'fontsize',14);


title('Solar Orbiter RPW');
