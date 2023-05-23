function quicklooks_7days(data,paths,Tint,logoPath)
% Given data in the struct 'data' (see solo.qli.quicklooks_main), generates
% plots and saves them in the paths specified in the struct 'paths' (see
% solo.qli.quicklooks_main). Tint should be a 7-day time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-08T00:00:00.00Z');

% Setup figure:
lwidth  = 1.0;
fsize   = 18;
legsize = 22;
h       = irf_plot(9,'newfigure');
fig     = gcf;
fig.Position = [1,1,1095,800];
colors       = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];


Units = irf_units;
Me    = Units.me;      % Electron mass [kg]
epso  = Units.eps0;    % Permitivitty of free space [Fm^-1]
mp    = Units.mp;      % Proton mass [km]
qe    = Units.e;       % Elementary charge [C]



%===================================
% Fill panel 1: B vector components
%===================================
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
    hold(h(1),'on');
    irf_plot(h(1),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',fsize);
irf_zoom(h(1),'y');



%======================
% Fill panel 2: abs(B)
%======================
if ~isempty(data.B)
    fci = qe*data.B.abs*10^-9/mp/(2*pi);
    irf_plot(h(2),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
ylabel(h(2),{'|B|';'(nT)'},'interpreter','tex','fontsize',fsize);
h(2).YScale='log';
h(2).YTick=[10,100];
%h(2).YLim=[0.1,200];



%=========================
% Fill panel 3: Densities
%=========================
hold(h(3),'on');
if ~isempty(data.Ne)
    irf_plot(h(3),data.Ne.tlim(Tint),'color',colors(1,:),'linewidth',lwidth);
end
if ~isempty(data.Npas)
    irf_plot(h(3),data.Npas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(3),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',fsize);
irf_legend(h(3),{'N_{e,RPW} ',' N_{i,PAS}'},[0.98 0.16],'Fontsize',legsize);
h(3).YScale='log';
h(3).YTick=[10,100];
%h(3).YLim=[0.8,200];



%===============================
% Fill panel 4: Ion temperature
%===============================
if ~isempty(data.Tpas)
    irf_plot(h(4),data.Tpas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(4),{'T_i';'(eV)'},'interpreter','tex','fontsize',fsize);
h(4).YScale='log';
h(4).YTick=[1,10,100];
h(4).YLim=[0.5,300];



%==============
% Fill panel 5
%==============
% y,z PAS velocities
if ~isempty(data.Vpas)
    irf_plot(h(5),data.Vpas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(5),'on');
    irf_plot(h(5),data.Vpas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(5),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(5),'y');
ylabel(h(5),{'v_{T,N}';'(km/s)'},'interpreter','tex','fontsize',fsize);



%==============
% Fill panel 6
%==============
hold(h(6),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(6),-data.Vrpw,'o','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(6),{'V_{RPW}','V_{PAS}'},[0.98 0.15],'Fontsize',legsize);
%h(6).YLim=[150,950];
ylabel(h(6),{'v_{R}';'(km/s)'},'interpreter','tex','fontsize',fsize);



%==============
% Fill panel 7
%==============
if ~isempty(data.E)
    irf_plot(h(7),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(7),'on');
    %irf_plot(h(7),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(7),{'','E_y'},[0.98 0.20],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',fsize);



%===================================
% Fill panel 8: Ion energy spectrum
%===================================
if ~isempty(data.ieflux)
    myFile=solo.db_list_files('solo_L2_swa-pas-eflux',Tint);
    iDEF   = struct('t',  data.ieflux.tlim(Tint).time.epochUnix);
    %for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
    for ii = 1:length(myFile)
        iEnergy = cdfread([myFile(ii).path '/' myFile(ii).name],'variables','Energy');
        iEnergy = iEnergy{1};
        iDEF.p = data.ieflux.data;
    end
    iDEF.f = repmat(iEnergy,1,numel(iDEF.t))';
    iDEF.p_label={'dEF','keV/','(cm^2 s sr keV)'};
    irf_spectrogram(h(8),iDEF,'log','donotfitcolorbarlabel');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    hold(h(8),'on');
    h8_clims = h(8).CLim;
    % Fix color axis
    h8_medp = mean(iDEF.p);
    h8_medp = min(h8_medp(h8_medp>0));
    if h8_medp > 0 && h8_medp > h8_clims(1) && log10(h8_medp)+2<(max(max(log10(iDEF.p))))
        caxis(h(8),[log10(h8_medp)+2 (max(max(log10(iDEF.p))))])
    end
    set(h(8), 'YScale', 'log');
    colormap(h(8),jet)
    ylabel(h(8),'[eV]')
end



%======================================
% Fill panel 9: E-field spectrum (TNR)
%======================================
if ~isempty(data.Etnr)
    % Electron plasma frequency
    myFile2=solo.db_list_files('solo_L2_rpw-tnr-surv-cdag',Tint);
    tp =[];pp=[];
    warning('off', 'fuzzy:general:warnDeprecation_Combine');
    TNR = [];
    %for iii = 1:round((myFile2(end).stop-myFile2(1).start)/3600/24)
    for iii = 1:length(myFile2)
        tt = [myFile2(iii).start myFile2(iii).stop];
        [TNRp] =  solo.read_TNR(tt);
        if isa(TNRp,'struct')
            TNR.t = combine(tp,TNRp.t);
            tp    = TNR.t;
            TNR.p = combine(pp,TNRp.p);
            pp    = TNR.p;

            % IMPLEMENTATION NOTE: Only read from TNRp from within this if
            % clause, since it might not be a struct if read from elsewhere,
            % even if it in principle means overwriting the value multiple
            % times as for TNRp.f and TNRp.p_label.
            TNR.f       = TNRp.f;
            TNR.p_label = TNRp.p_label;
        end
    end
    if isstruct(TNR)
        % TNR.f       = TNRp.f;
        % TNR.p_label = TNRp.p_label;
        sz_tnr = size(TNR.p);
        if sz_tnr(1) == length(TNR.t) && sz_tnr(2) == length(TNR.f)
            irf_spectrogram(h(9),TNR,'log','donotfitcolorbarlabel')
            hold(h(9),'on');
        end
        if ~isempty(data.Ne)
            wpe_sc = (sqrt(((data.Ne.tlim(Tint)*1000000)*qe^2)/(Me*epso)));
            fpe_sc = (wpe_sc/2/pi)/1000;
            fpe_sc.units = 'kHz';
            fpe_sc.name  = 'f [kHz]';
            irf_plot(h(9),fpe_sc,'r','linewidth',lwidth);
        end
        text(h(9),0.01,0.3,'f_{pe,RPW}','units','normalized','fontsize',18,'Color','r');
        %set(h(9), 'YScale', 'log');
        colormap(h(9),jet)
        %ylabel(h(9),'f [kHz]')
        set(h(9),'ColorScale','log')
        %caxis([.01 10]*10^-12)
        yticks(h(9),[10^1 10^2]);
    end
end



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:9));
irf_zoom(h(1:9),'x',Tint);
irf_zoom(h(1),'y');
irf_zoom(h(5:9),'y');

h(2).YLabel.Position=[1.05,0.5,0];
%yyaxis(h(2),'left');
h(2).YLabel.Units='normalized';
h(2).YLabel.Position=h(3).YLabel.Position;
h(9).XLabel.Visible = 'off';



% if ~isempty(data.solopos.tlim(Tint))
%     teststr = ['SolO: ', ...
%         [' R=', sprintf('%.2f',data.solopos.tlim(Tint   ).data(1,1)/AU_KM),'Au, '],...
%         [' EcLat=',sprintf('%d',round(data.solopos.tlim(Tint   ).data(1,3)*180/pi)),'\circ, '],...
%         [' EcLon=',sprintf('%d',round(data.solopos.tlim(Tint   ).data(1,2)*180/pi)),'\circ']];
%     text1=text(h(9),-0.11,-0.575,teststr,'units','normalized','fontsize',18);
% else
%     teststr = char();
%     text1   = text(h(9),-0.11,-0.575,teststr,'units','normalized','fontsize',18);
% end
[soloStr, earthStr] = solo.qli.context_info_strings(data.solopos, data.earthpos, Tint);
text(h(9), -0.11, -0.575, soloStr, 'units', 'normalized', 'fontsize', 18);

% Add Earth longitude as text.
% if ~isempty(data.earthpos)
%     teststr =['Earth: EcLon=',sprintf('%d',round(data.earthpos.data(1,2)*180/pi)),'\circ'];
%     text2=text(h(9),-0.11,-0.925,teststr,'units','normalized','fontsize',18);
% else
%     teststr=char();
%     text2=text(h(9),-0.11,-0.925,teststr,'units','normalized','fontsize',18);
% end
text(h(9), -0.11, -0.925, earthStr, 'units', 'normalized', 'fontsize', 18);



xtickangle(h(9),0)
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
currdate = char(datetime("now","Format","uuuu-MM-dd"));
infostr = ['Swedish Institute of Space Physics, Uppsala (IRFU), ',currdate];
infostr2 = '. Data available at http://soar.esac.esa.int/';
text(h(1),0,1.2,[infostr,infostr2],'Units','normalized')

% Remove overlapping tics.
solo.qli.ensure_axes_data_tick_margins(h)

%yyaxis(h(2),'left');
oldlims2 = h(2).YLim;
oldticks2 = h(2).YTick;
h(2).YScale='log';
h(2).YTick=[1,10,100];
h(2).YLim=[0.8,200];

% yyaxis(h(2),'right');
% oldlims2_r=h(2).YLim;
% oldticks2_r = h(2).YTick;
% h(2).YScale='log';
% h(2).YTick=[1,10,100];
%h(2).YLim=[0.1,200];

oldlims5 = h(5).YLim;
oldticks5 = h(5).YTick;
h(5).YScale='log';
h(5).YTick=[1,10,100];
%h(5).YLim=[0.5,300];

c_eval('h(?).FontSize=18;',1:9);


irf_plot_axis_align(h(1:9));
irf_zoom(h(1:9),'x',Tint);
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
%=====================
% Save figure to file
%=====================
print('-dpng',path1);
% TODO-NI: Why are there any commands (except close()) after this?

% h(2).YScale='lin';
% h(2).YTick=oldticks2_r;
% h(2).YLim=oldlims2_r;
%yyaxis(h(2),'left');
h(2).YScale='lin';
h(2).YLim=oldlims2;
h(2).YTick=oldticks2;

h(5).YScale='lin';
h(5).YLim=oldlims5;
h(5).YTick=oldticks5;

close(fig);

end
