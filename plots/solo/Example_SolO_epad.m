%% Read data
tint = irf.tint('2023-03-20T00:00:00.000Z/2023-03-23T00:00:00.000Z');
b = solo.get_data('b_rtn_norm',tint); % read magnetic field
v = solo.get_data('vi_rtn',tint); % read velocity
n = solo.get_data('ni',tint); % read density
ePAD = solo.get_data('epad_10sec',tint); % read electron pad data

epad_sepc = ePAD.ebin_32; % take data from index 32 (see ePAD.enegies for the energy bin values)

% add nans at the edges of datagaos to stop linear interpolation
gap_ind = find(diff(epad_sepc.t)>2*mode(diff(epad_sepc.t)));
epad_sepc.p(gap_ind,:) = nan;
epad_sepc.p(gap_ind+1,:) = nan;
gap_ind = find(diff(n.time.epochUnix)>2*mode(diff(n.time.epochUnix)));
n.data(gap_ind,:) = nan;
n.data(gap_ind+1,:) = nan;
v.data(gap_ind,:) = nan;
v.data(gap_ind+1,:) = nan;


%% plot data
h1 = irf_figure(5);
irf_plot(h1(1),b.abs)
irf_plot(h1(2),b)
irf_plot(h1(3),v)
irf_plot(h1(4),n)
irf_spectrogram(h1(5),epad_sepc,'log','donotshowcolorbar')
hcb2 = colorbar(h1(5));
ylabel(hcb2,{'log10 PSD' ;'[s^3 km^{-6}]'},'fontsize',15)
clim(h1(5),[2 4.5])
irf_plot_axis_align(h1);
colormap jet
set(gcf,'Position',[1800 300 1100 900])
set(h1, 'TickDir', 'out','linewidth',2,'FontSize',15);
set(h1(5) ,'Layer', 'Top')
clim([2.5 4])
ylabel(h1(1),'|B| [nT]','FontSize',15)
ylabel(h1(2),'Brtn [nT]','FontSize',15)
ylabel(h1(3),'Vrtn [km/s]','FontSize',15)
ylabel(h1(4),'ni [cc]','FontSize',15)
ylabel(h1(5),'Angle [deg]','FontSize',15)
grid(h1,'off')
irf_zoom(h1,'x',tint)

ylim(h1(1),[0 80])
set(h1(1),'YTick',[0:10:80]);

ylim(h1(2),[-60 60])
set(h1(2),'YTick',[-40:20:40]);

ylim(h1(3),[-300 1000])
set(h1(3),'YTick',[-200:200:800]);

ylim(h1(4),[0 220])
set(h1(4),'YTick',[0:25:200]);

ylim(h1(5),[0 180])
set(h1(5),'YTick',[0:30:150]);


irf_legend(h1(2),{'Br';'Bt';'Bn'},[1.02 0.98],'fontsize',15,'FontName','monospaced')
irf_legend(h1(3),{'Vr';'Vt';'Vn'},[1.02 0.98],'fontsize',15,'FontName','monospaced')
irf_legend(h1(5),{['Electrons: ' num2str(round(ePAD.energies(32))) ' eV bin']},[0.99 0.1],'fontsize',15,'FontName','monospaced','backgroundcolor','w')

irf_legend(h1(1),{'(a)'},[0.99 0.95],'fontsize',15,'FontName','monospaced','backgroundcolor','w')
irf_legend(h1(2),{'(b)'},[0.99 0.95],'fontsize',15,'FontName','monospaced','backgroundcolor','w')
irf_legend(h1(3),{'(c)'},[0.99 0.95],'fontsize',15,'FontName','monospaced','backgroundcolor','w')
irf_legend(h1(4),{'(d)'},[0.99 0.95],'fontsize',15,'FontName','monospaced','backgroundcolor','w')
irf_legend(h1(5),{'(e)'},[0.99 0.95],'fontsize',15,'FontName','monospaced','backgroundcolor','w')

colormap(h1(5),irf_colormap('waterfall'))
clim(h1(5),[3 4])
hcb2.LineWidth = 2;