% plot figure showing typical temporal and spatial scales and doppler
% shifts for solar orbiter plasma conditions

NaturalConstants

f_range=[0.1 1e6]; % frequencies in Hz
l_range=[1 1e6]; % wave length in m 
l_inv_range=fliplr(1./l_range); % inverse wave length
v_kms=[100 300 1000]; % velocity lines in km/s
v=v_kms*1e3; % to get in m (SI units)

figure(1);clf;
h=subplot(1,1,1);
set(gca,'xscale','log','xlim',l_inv_range);
set(gca,'yscale','log','ylim',f_range);
color_order=[0 0 0;1 0 0;0 0 1];
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
    ht=text(l_inv_range(1)*2,f_range(2)*0.95*0.5^jj,ttt);
    set(ht,'FontSize',14,'VerticalAlignment','top','Color',color_order(jj,:));
end

patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [f_range(1) f_range(1) 300 300],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.12, 10,'DC');set(ht,'fontsize',14);
patch(l_inv_range(2)*[0.33 0.99 0.99 0.33], [100 100 25000 25000],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.4, 1000,'LF');set(ht,'fontsize',14);
patch(l_inv_range(2)*[0.1 0.3 0.3 0.1], [1e4 1e4 1e6 1e6],[1 0.5 0.5]) 
ht=text(l_inv_range(2)*0.12, 1e5,'HF');set(ht,'fontsize',14);

title('Solar Orbiter RPW');
