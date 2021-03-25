function quicklooks_7days(data,paths,Tint)
% Given data in the struct 'data' (see solo.quicklook_main), generates
% plots and saves in the paths specified in the struct 'paths' (see
% solo.quicklook_main). Tint should be a 7-day time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-08T00:00:00.00Z');

% Setup figure:
lwidth=1.0;
fsize=18;
legsize=22;
h=irf_plot(8,'newfigure');
fig=gcf;
fig.Position=[1,1,1095,800];
colors = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),'B_{RTN} (nT)','interpreter','tex','fontsize',fsize);

if ~isempty(data.B)
    irf_plot(h(2),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
ylabel(h(2),'|B| (nT)','interpreter','tex','fontsize',fsize);

% Densities
hold(h(3),'on');
if ~isempty(data.Ne)
    irf_plot(h(3),data.Ne.tlim(Tint),'color',colors(1,:),'linewidth',lwidth);
else
end
if ~isempty(data.Npas)
    irf_plot(h(3),data.Npas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(3),'N (cm^{-3})','interpreter','tex','fontsize',fsize);
irf_legend(h(3),{'N_{e,RPW} ',' N_{i,PAS}'},[0.98 0.16],'Fontsize',legsize);
irf_zoom(h(3),'y');

if ~isempty(data.Tpas)
    irf_plot(h(4),data.Tpas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(4),'T_i (eV)','interpreter','tex','fontsize',fsize);
irf_zoom(h(4),'y');

% y,z PAS velocities
if ~isempty(data.Vpas)
    irf_plot(h(5),data.Vpas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(5),'on');
    irf_plot(h(5),data.Vpas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(5),{'','v_{y}','v_{z}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(5),'y');
ylabel(h(5),'v_{yz} (km/s)','interpreter','tex','fontsize',fsize);

hold(h(6),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(6),data.Vrpw,'o-','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(6),{'V_{RPW}','V_{PAS}'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(6),'y');
ylabel(h(6),'v_{x} (km/s)','interpreter','tex','fontsize',fsize);

if ~isempty(data.E)
    irf_plot(h(7),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(7),'on');
    irf_plot(h(7),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(7),{'','E_y','E_z'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),'E (mV/m)','interpreter','tex','fontsize',fsize);

Au=149597871;
irf_plot(h(8),data.solopos.x/Au);

ylabel(h(8),'R_{sun} (Au)','interpreter','tex','fontsize',fsize);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
irf_zoom(h(1:8),'y');

% Plot complete, print figure.
fig=gcf;
fig.PaperPositionMode='auto';

filesmth = Tint(1);
filesmth = filesmth.utc;
filestr1 = filesmth(1:13);
filestr1([5,8])=[];

filesmth = Tint(end);
filesmth = filesmth.utc;
filestr2 = filesmth(1:13);
filestr2([5,8])=[];
path1=fullfile(paths.path_1w,[filestr1,'_',filestr2,'.png']);
print('-dpng',path1);

close(fig);
