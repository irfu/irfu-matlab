%% Load the data and make the PDist object
Tint = irf.tint('2021-08-22T00:00:00.00Z/2021-08-23T00:00:00.00Z'); % Tint for the entire day

PDout = solo.get_data('pas_vdf',Tint); % read the SWA-PAS VDF
PDout = PDout.(irf_time(Tint.start,'epochtt>utc_Tyyyymmdd')); % Take the PDist for this day
SCpot = irf.ts_scalar(PDout.time,zeros(size(PDout.time)));

%B = solo.get_data('B_rtn_norm',Tint);
%% Distributions
Tint_vdf = irf.tint('2021-08-22T02:19:30.00Z/2021-08-22T02:21:30.00Z'); % Time you want the VDF

PDout2 = PDout.tlim(Tint_vdf);
SCpot = SCpot.tlim(Tint_vdf);
Ng = 100; nMC = 1e3; vg = linspace(-15e2,15e2,Ng);
nAverage = 10;
indT = 1:nAverage; % Plots average of first nAverage distributions in Tint_vdf interval

%Bav = B.tlim(Tint_vdf); Bav = median(double(Bav.data)); nB = irf_norm(Bav);

f2Dxy = PDout2(indT).reduce('2D',[1 0 0],[0 1 0],'base','cart','vg',vg,'nMC',nMC,'SCpot',SCpot,'vint',[-1500 1500]);
f2Dxz = PDout2(indT).reduce('2D',[1 0 0],[0 0 1],'base','cart','vg',vg,'nMC',nMC,'SCpot',SCpot,'vint',[-1500 1500]);
f2Dyz = PDout2(indT).reduce('2D',[0 1 0],[0 0 1],'base','cart','vg',vg,'nMC',nMC,'SCpot',SCpot,'vint',[-1500 1500]);

%% Plot distributions
fn=figure;
xwidth = 0.25;
ywidth = 0.75;

set(fn,'Position',[10 10 1100 300])
h(1)=axes('position',[0.06 0.18 xwidth ywidth]); % [x y dx dy]
h(2)=axes('position',[0.39 0.18 xwidth ywidth]);
h(3)=axes('position',[0.72 0.18 xwidth ywidth]);
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',1);

f2Dxy.plot_plane(h(1),'docolorbar',0);
caxis(h(1),[-7   -1]);
hcb = colorbar(h(1));
ylabel(hcb,'log_{10} f_i (s^2 m^{-5})','interpreter','tex','fontsize',12)
grid(h(1),'off')
axis(h(1),'equal')
colormap(h(1),'jet')
irf_legend(h(1),'(a)',[0.98 0.98],'color','k','fontsize',12)
axis(h(1),[-600 -100 -250 250])
xlabel(h(1),'V_{x} (km s^{-1})','interpreter','tex','fontsize',12)
ylabel(h(1),'V_{y} (km s^{-1})','interpreter','tex','fontsize',12)
if nAverage>1, title(h(1),f2Dxy.time([1 nAverage]).toUtc)
else, title(h(1),f2Dxy(1).time.toUtc)
end

f2Dxz.plot_plane(h(2),'docolorbar',0);
caxis(h(2),[-7   -1]);
hcb = colorbar(h(2));
ylabel(hcb,'log_{10} f_i (s^2 m^{-5})','interpreter','tex','fontsize',12)
grid(h(2),'off')
axis(h(2),'equal')
colormap(h(2),'jet')
irf_legend(h(2),'(b)',[0.98 0.98],'color','k','fontsize',12)
axis(h(2),[-600 -100 -250 250])
xlabel(h(2),'V_{x} (km s^{-1})','interpreter','tex','fontsize',12)
ylabel(h(2),'V_{z} (km s^{-1})','interpreter','tex','fontsize',12)
if nAverage>1, title(h(2),f2Dxz.time([1 nAverage]).toUtc)
else, title(h(2),f2Dxz(1).time.toUtc)
end

f2Dyz.plot_plane(h(3),'docolorbar',0);
caxis(h(3),[-7   -1]);
hcb = colorbar(h(3));
ylabel(hcb,'log_{10} f_i (s^2 m^{-5})','interpreter','tex','fontsize',12)
grid(h(3),'off')
axis(h(3),'equal')
colormap(h(3),'jet')
irf_legend(h(3),'(c)',[0.98 0.98],'color','k','fontsize',12)
axis(h(3),[-250 250 -250 250])
xlabel(h(3),'V_{y} (km s^{-1})','interpreter','tex','fontsize',12)
ylabel(h(3),'V_{z} (km s^{-1})','interpreter','tex','fontsize',12)
if nAverage>1, title(h(3),f2Dyz.time([1 nAverage]).toUtc)
else, title(h(3),f2Dyz(1).time.toUtc)
end

set(h(1:3),'fontsize',12)
set(gcf,'color','w');