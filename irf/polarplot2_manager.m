function status=polarplot2_manager(arg)
persistent flag_new_figure h2 dgud;
global hf
global POLARPLOT_RESULTS   % POLARPLOT_RESULTS={t0, T-t0,freq,2*n*Imean1/sampl,coh1,phase1,2*n*Imean2/sampl,coh2,phase2,threshold};
t0=POLARPLOT_RESULTS{1};
t=POLARPLOT_RESULTS{2};
freq=POLARPLOT_RESULTS{3};
S1=POLARPLOT_RESULTS{4};
Coh1=POLARPLOT_RESULTS{5};
Phase1=POLARPLOT_RESULTS{6};
S2=POLARPLOT_RESULTS{7};
Coh2=POLARPLOT_RESULTS{8};
Phase2=POLARPLOT_RESULTS{9};
threshold=POLARPLOT_RESULTS{10};

p = get(gca, 'currentpoint');
phase_lim=get(gca,'Clim');

ind=find(p(1)>t, 1, 'last');
ss1=S1(:,ind);
ph1=Phase1(:,ind);
ch1=Coh1(:,ind);
ss2=S2(:,ind);
ph2=Phase2(:,ind);
ch2=Coh2(:,ind);

if isempty(flag_new_figure)
  freq_lim=get(gca,'Ylim');
  hf=figure;flag_new_figure=1;clf;
else
  freq_lim=[str2num(get(dgud.fplminh,'string')) str2num(get(dgud.fplmaxh,'string')) ];
  figure(hf);
end

if flag_new_figure==1, h2(1)=irf_subplot(4,1,-1); end
axes(h2(1));
semilogy(freq,ss1,'k',freq,ss2,'b');set(gca,'xlim',freq_lim);grid on;
ylabel('Power [(mV/m)^2/Hz]');
title([datestr(datenum(fromepoch(p(1)+t0))) '.' num2str(mod(p(1)+t0,1)*100,2) ' Spectra']);

if flag_new_figure==1, h2(2)=irf_subplot(4,1,-2);end
axes(h2(2));
plot(freq,ch1,'k.',freq,ch2,'b.');axis([freq_lim 0 1]);axis manual;
set(gca,'ytick',[0 threshold 1]);grid on;
ylabel('Coherence');
ind_coh_good1=find(ch1>threshold);
ind_coh_bad1=find(ch1<=threshold);
ind_coh_good2=find(ch2>threshold);
ind_coh_bad2=find(ch2<=threshold);

if flag_new_figure==1, h2(3)=irf_subplot(4,1,-3);end
axes(h2(3));cla;
plot(freq(ind_coh_good1),ph1(ind_coh_good1),'k.',freq(ind_coh_bad1),ph1(ind_coh_bad1),'kx');
hold on;grid on;
plot(freq(ind_coh_good2),ph2(ind_coh_good2),'b.',freq(ind_coh_bad2),ph2(ind_coh_bad2),'bx');
axis([freq_lim phase_lim]);%set(gca,'ytick',-180:45:180);
hl1=line([0 0],[0 0]);set(hl1,'color','k');
hl2=line([0 0],[0 0]);set(hl2,'color','b');
ylabel('Phase');

if flag_new_figure==1
  h2(4)=subplot(4,1,4);axis off;

  dgud.h2=h2;
  ysc=1.5;% scaling for y separation between text information on screen
  dgud.fmint=uicontrol('style', 'text', 'string', 'fmin  [Hz]','units','centimeters','position', [2 ysc*1 3 ysc*.5]);
  dgud.fminh = uicontrol('style', 'edit', 'string','0', 'units','centimeters','position', [5 ysc*1 3 ysc*.5],'callback', 'polarplot2_manager_v(''recalculate'')');
  dgud.fmax=uicontrol('style', 'text', 'string', 'fmax [Hz]','units','centimeters','position', [2 ysc*1.5 3 ysc*.5]);
  dgud.fmaxh = uicontrol('style', 'edit', 'string', '0', 'units','centimeters','position', [5 ysc*1.5 3 ysc*.5], 'callback','polarplot2_manager_v(''recalculate'')');
  dgud.fplmin=uicontrol('style', 'text', 'string', 'fplot min [Hz]','units','centimeters','position', [2 ysc*2 3 ysc*.5]);
  dgud.fplminh = uicontrol('style', 'edit', 'string', num2str(freq_lim(1)),'units','centimeters','position', [5 ysc*2 3 ysc*.5],'callback', 'polarplot2_manager_v(''fpl'')');
  dgud.fplmax=uicontrol('style', 'text', 'string', 'fplot max [Hz]', 'units','centimeters','position', [2 ysc*2.5 3 ysc*.5]);
  dgud.fplmaxh = uicontrol('style', 'edit', 'string', num2str(freq_lim(2)), 'units','centimeters','position', [5 ysc*2.5 3 ysc*.5],'callback', 'polarplot2_manager_v(''fpl'')');

  dgud.vfitzero=uicontrol('style', 'checkbox', 'string', 'go through zero', 'units','centimeters','position', [9 ysc*2.5 3 ysc*.5],'callback', 'polarplot2_manager_v(''recalculate'')');
  dgud.sampling_distance=uicontrol('style', 'text', 'string', 'distance in phys.un.', 'units','centimeters','position', [9 ysc*2 5 ysc*.5]);
  dgud.sampling_distance_h=uicontrol('style', 'edit', 'string', '0.044','units','centimeters','position', [14 ysc*2 3 ysc*.5],'callback', 'polarplot2_manager_v(''recalculate'')');
  dgud.vphase1=uicontrol('style', 'text', 'string', '','units','centimeters','position', [9 ysc*1.5 8 ysc*.5]);
  dgud.vphase2=uicontrol('style', 'text', 'string', '','units','centimeters','position', [9 ysc*1   8 ysc*.5]);
end

dgud.fitline1=hl1;
dgud.fitline2=hl2;

dgud.ph1=ph1;dgud.ph2=ph2;dgud.freq=freq;
dgud.ind_coh_good1=ind_coh_good1;dgud.ind_coh_good2=ind_coh_good2;

set(gcf,'userdata',dgud);

if flag_new_figure==0,polarplot2_manager_v('recalculate');end

flag_new_figure=0;


