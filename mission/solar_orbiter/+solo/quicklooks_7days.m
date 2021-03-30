function quicklooks_7days(data,paths,Tint,logoPath)
% Given data in the struct 'data' (see solo.quicklook_main), generates
% plots and saves in the paths specified in the struct 'paths' (see
% solo.quicklook_main). Tint should be a 7-day time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-08T00:00:00.00Z');

% Setup figure:
lwidth=1.0;
fsize=18;
legsize=22;
h=irf_plot(7,'newfigure');
fig=gcf;
fig.Position=[1,1,1095,800];
colors = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
    hold(h(1),'on');
    irf_plot(h(1),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',fsize);

if ~isempty(data.B)
    irf_plot(h(2),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
ylabel(h(2),{'|B|';'(nT)'},'interpreter','tex','fontsize',fsize);
h(2).YScale='log';
h(2).YTick=[1,10,100];
h(2).YLim=[0.1,200];

% Densities
hold(h(3),'on');
if ~isempty(data.Ne)
    irf_plot(h(3),data.Ne.tlim(Tint),'color',colors(1,:),'linewidth',lwidth);
else
end
if ~isempty(data.Npas)
    irf_plot(h(3),data.Npas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(3),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',fsize);
irf_legend(h(3),{'N_{e,RPW} ',' N_{i,PAS}'},[0.98 0.16],'Fontsize',legsize);
h(3).YScale='log';
h(3).YTick=[1,10,100];
h(3).YLim=[0.8,200];

if ~isempty(data.Tpas)
    irf_plot(h(4),data.Tpas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(4),{'T_i';'(eV)'},'interpreter','tex','fontsize',fsize);
h(4).YScale='log';
h(4).YTick=[1,10,100];
h(4).YLim=[0.5,300];


% y,z PAS velocities
if ~isempty(data.Vpas)
    irf_plot(h(5),data.Vpas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(5),'on');
    irf_plot(h(5),data.Vpas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(5),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(5),'y');
ylabel(h(5),{'v_{T,N}';'(km/s)'},'interpreter','tex','fontsize',fsize);

hold(h(6),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(6),-data.Vrpw,'o','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(6),{'V_{RPW}','V_{PAS}'},[0.98 0.15],'Fontsize',legsize);
h(6).YLim=[150,950];
ylabel(h(6),{'v_{R}';'(km/s)'},'interpreter','tex','fontsize',fsize);

if ~isempty(data.E)
    irf_plot(h(7),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(7),'on');
    irf_plot(h(7),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(7),{'','E_y','E_z'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',fsize);

Au=149597871; %Astronomical unit.

if ~isempty(data.solopos.tlim(Tint))
    teststr = ['SolO: ',[' R=',sprintf('%.2f',data.solopos.tlim(Tint).data(1,1)/Au),'Au, '],...
        [' EcLat=',sprintf('%d',round(data.solopos.tlim(Tint).data(1,3)*180/pi)),'°, '],...
        [' EcLon=',sprintf('%d',round(data.solopos.tlim(Tint).data(1,2)*180/pi)),'°']];
    text1=text(h(7),-0.11,-0.4,teststr,'units','normalized','fontsize',18);


else
    teststr=char();
    text1=text(h(7),-0.11,-0.4,teststr,'units','normalized','fontsize',18);

end

% Add Earth longitude as text.
if ~isempty(data.earthpos)
    teststr =['Earth: EcLon=',sprintf('%d',round(data.earthpos.data(1,2)*180/pi)),'°'];
    text2=text(h(7),-0.11,-0.65,teststr,'units','normalized','fontsize',18);
else
    teststr=char();
    text2=text(h(7),-0.11,-0.65,teststr,'units','normalized','fontsize',18);
end

% Add plot information and IRF logo
logopos = h(1).Position;
logopos(1)=logopos(1)+logopos(3)+0.01;
logopos(2)=logopos(2)+0.06;
logopos(3)=0.05;
logopos(4)=logopos(3)*1095/800;
ha2=axes('position',logopos);

if ~isempty(logoPath)
    [x, map]=imread(logoPath);
    image(x)
end
% colormap (map)
set(ha2,'handlevisibility','off','visible','off')
tempdate=datestr(date,2);
currdate=['20',tempdate(7:8),'-',tempdate(1:2),'-',tempdate(4:5)];
infostr = ['Swedish Institute of Space Physics, Uppsala (IRFU), ',currdate];
infostr2 = '. Data available at http://soar.esac.esa.int/';
text(h(1),0,1.2,[infostr,infostr2],'Units','normalized')

c_eval('h(?).FontSize=18;',1:7);


irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',Tint);
% irf_zoom(h(1:7),'y');

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
