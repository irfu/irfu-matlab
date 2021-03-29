function quicklooks_24_6_2_h(data,paths,Tint)
% Given data in the struct 'data' (see solo.quicklook_main), generates
% plots and saves in the paths specified in the struct 'paths' (see
% solo.quicklook_main). Computes spectrum of B, so takes a while to run.
% Tint should be a 24hour time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-02T00:00:00.00Z');


% Setup figure:
lwidth=1.0;
fsize=18;
legsize=22;
h=irf_plot(8,'newfigure');
fig=gcf;
fig.Position=[1,1,1095,800];
colors = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];
mp=1.6726e-27; %proton mass [km]
qe=1.6022e-19; %elementary charge [C]
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
    hold(h(1),'on');
    irf_plot(h(1),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',fsize);

%%
hold(h(2),'on');
if ~isempty(data.Ne)
    irf_plot(h(2),data.Ne.tlim(Tint),'-','color',colors(1,:),'linewidth',lwidth);
end
if ~isempty(data.Npas)
    irf_plot(h(2),data.Npas.tlim(Tint),'-','color',colors(2,:),'linewidth',lwidth);
end
ylb2=ylabel(h(2),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',fsize);
h(2).ColorOrder=colors;
irf_legend(h(2),{'N_{e,RPW}','N_{i,PAS}','|B|'},[0.98 0.16],'Fontsize',legsize);
irf_zoom(h(2),'y');

yyaxis(h(2),'right');
if ~isempty(data.B)
    fci = qe*data.B.abs*10^-9/mp/(2*pi);
    irf_plot(h(2),data.B.abs.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
    h(2).YLim=[floor(min(data.B.abs.data)),ceil(max(data.B.abs.data))];
end
ylabel(h(2),{'|B|';'(nT)'},'interpreter','tex','fontsize',fsize);
h(2).YColor=[1,0,0];


%%

%Spectra h(3), h(4)
if ~isempty(data.B)
    bb = data.B;
    if median(diff((bb.time.epochUnix))) < 0.1250*0.95
        fMag = 128; fMax = 7;
    else
        fMag = 8; fMax = 3;
    end
    b0 = bb.filt(0, 0.01,fMag, 5);
    ebsp = irf_ebsp([],bb,[],b0,[],[0.05 fMax],'fullB=dB', 'polarization', 'fac');
    
    frequency = ebsp.f;
    time = ebsp.t;
    Bsum = ebsp.bb_xxyyzzss(:,:,4);
    ellipticity = ebsp.ellipticity;
    dop = ebsp.dop;
    
    % Remove points with very low degree of polarization
    dopthresh = 0.7;
    removepts = find(dop < dopthresh);
    ellipticity(removepts) = NaN;
    
    % Remove "lonely" pixels
    msk = ellipticity;
    msk(~isnan(msk)) = 1;
    msk(isnan(msk)) = 0;
    msk_denoise = bwareaopen(msk,8);
    ellipticity(msk_denoise==0) = NaN;
    
    % Plot
    specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=Bsum;
    specrec.f_label='';
    specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
    irf_spectrogram(h(3),specrec,'log','donotfitcolorbarlabel');
    set(h(3),'yscale','log');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    % caxis(h(3),[-8 -1])
    hold(h(3),'on');
    irf_plot(h(3),fci,'k','linewidth',lwidth);
    text(h(3),0.01,0.3,'f_{ci}','units','normalized','fontsize',18);
    ylabel(h(3),{'f';'(Hz)'},'fontsize',fsize);
    colormap(h(3),'jet');
    
    
    specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=ellipticity;
    specrec.f_label='';
    specrec.p_label={'Ellipticity','DOP>0.7'};
    irf_spectrogram(h(4),specrec,'log','donotfitcolorbarlabel');
    set(h(4),'yscale','log');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    caxis(h(4),[-1 1])
    hold(h(4),'on');
    irf_plot(h(4),fci,'k','linewidth',lwidth);
    ylabel(h(4),{'f';'(Hz)'},'fontsize',fsize);
    text(h(4),0.01,0.3,'f_{ci}','units','normalized','fontsize',18);
    
    crr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
    cgg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
    cbb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
    bgrcmap = [crr' cgg' cbb'];
    colormap(h(4),bgrcmap);
end


% Ion temperature
if ~isempty(data.Tpas)
    irf_plot(h(5),data.Tpas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(5),{'T_i';'(eV)'},'interpreter','tex','fontsize',fsize);
irf_zoom(h(5),'y');

% y,z PAS velocities
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(6),'on');
    irf_plot(h(6),data.Vpas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(6),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(6),'y');
ylabel(h(6),{'V_{T,N}';'(km/s)'},'interpreter','tex','fontsize',fsize);

hold(h(7),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(7),-data.Vrpw,'o-','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(7),data.Vpas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(7),{'V_{RPW}','V_{PAS}'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),{'V_R';'(km/s)'},'interpreter','tex','fontsize',fsize);

if ~isempty(data.E)
    irf_plot(h(8),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(8),'on');
    irf_plot(h(8),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(8),{'','E_y','E_z'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(8),'y');
ylabel(h(8),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',fsize);



irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
irf_zoom(h(1),'y');
irf_zoom(h(5:8),'y');
h(2).YLabel.Position=[1.05,0.5,0];
yyaxis(h(2),'left');
h(2).YLabel.Units='normalized';
h(2).YLabel.Position=h(3).YLabel.Position;
% Add spacecraft position as text.
Au=149597871; %km
if ~isempty(data.solopos.tlim(Tint))
    teststr = ['SolO: ',[sprintf('%.2f',data.solopos.tlim(Tint).data(1,1)/Au),'Au, '],...
        [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint).data(1,3)*180/pi)),'°, '],...
        [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint).data(1,2)*180/pi)),'°']];
    text1=text(h(8),-0.11,-0.5,teststr,'units','normalized','fontsize',18);
else
    teststr=char();
    text1=text(h(8),-0.11,-0.5,teststr,'units','normalized','fontsize',18);
end

% Add Earth longitude as text.
if ~isempty(data.earthpos)
    teststr =['Earth: EcLon ',sprintf('%d',round(data.earthpos(1,2)*180/pi)),'°'];
    text2=text(h(8),-0.11,-0.75,teststr,'units','normalized','fontsize',18);
else
    teststr=char();
    text2=text(h(8),0.9,-0.5,teststr,'units','normalized','fontsize',18);
end


% Add plot information and IRF logo
logopos = h(1).Position;
logopos(1)=logopos(1)+logopos(3)+0.04;
logopos(2)=logopos(2)+0.06;
logopos(3)=0.05;
logopos(4)=logopos(3)*1095/800;
ha2=axes('position',logopos);

[x, map]=imread('irf_logo.png');
image(x)
% colormap (map)
set(ha2,'handlevisibility','off','visible','off')
tempdate=datestr(date,2);
currdate=['20',tempdate(7:8),'-',tempdate(1:2),'-',tempdate(4:5)];
infostr = ['Swedish Institute of Space Physics, Uppsala (IRFU), ',currdate];
infostr2 = '. Data available at http://soar.esac.esa.int/';
text(h(1),0,1.2,[infostr,infostr2],'Units','normalized')

% Fix YTicks
for iax=1:8
    cax=h(iax);
    mintick = min(cax.YTick);
    maxtick = max(cax.YTick);
    minlim = cax.YLim(1);
    maxlim = cax.YLim(2);
    
    if maxtick>0
        if maxlim<1.1*maxtick
            newmax = 1.1*maxtick;
        else
            newmax = maxlim;
        end  
    else
        if abs(maxlim)>0.9*abs(maxtick)
            newmax = 0.9*maxtick;
        else
            newmax = maxlim;
        end       
    end
    
    if mintick>0  
        if minlim>0.9*mintick
            newmin = 0.9*mintick;
        else
            newmin = minlim;
        end
    else
        if abs(minlim)<1.1*abs(mintick)
            newmin=1.1*mintick;
        else
            newmin=minlim;
        end
    end
    cax.YLim=[newmin,newmax];
end

yyaxis(h(2),'left');
oldlims2 = h(2).YLim;
oldticks2 = h(2).YTick;
h(2).YScale='log';
h(2).YTick=[1,10,100];
h(2).YLim=[0.8,200];

yyaxis(h(2),'right');
oldlims2_r=h(2).YLim;
oldticks2_r = h(2).YTick;
h(2).YScale='log';
h(2).YTick=[1,10,100];
h(2).YLim=[0.1,200];

oldlims5 = h(5).YLim;
oldticks5 = h(5).YTick;
h(5).YScale='log';
h(5).YTick=[1,10,100];
h(5).YLim=[0.5,300];

% Plot complete. Print in 24h, 6h and 2h intervals.
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
path1=fullfile(paths.path_24h,[filestr1,'_',filestr2,'.png']);
print('-dpng',path1);

h(2).YScale='lin';
h(2).YTick=oldticks2_r;
h(2).YLim=oldlims2_r;
yyaxis(h(2),'left');
h(2).YScale='lin';
h(2).YLim=oldlims2;
h(2).YTick=oldticks2;



h(5).YScale='lin';
h(5).YLim=oldlims5;
h(5).YTick=oldticks5;
% Print 6h figures
for i6h = 1:4
    
    %Zoom in to 6h interval and save plot.
    Tint_6h = Tint(1)+[60*60*6*(i6h-1),60*60*6*(i6h)];
    irf_zoom(h(1:8),'x',Tint_6h);
    irf_zoom(h(1),'y');
    %Zoom on N/|B| plot..
    Neflag = ~isempty(data.Ne.tlim(Tint_6h));
    Npasflag = ~isempty(data.Npas) && ~isempty(data.Npas.tlim(Tint_6h));
    if Neflag && Npasflag
        yyaxis(h(2),'left');
        h(2).YLim=[min(floor([min(data.Npas.tlim(Tint_6h).data),min(data.Ne.tlim(Tint_6h).data)])),...
            max(ceil([max(data.Npas.tlim(Tint_6h).data),max(data.Ne.tlim(Tint_6h).data)]))];
    elseif Neflag
        yyaxis(h(2),'left');
        h(2).YLim=[floor(min(data.Ne.tlim(Tint_6h).data)),ceil(max(data.Ne.tlim(Tint_6h).data))];
    elseif Npasflag
        yyaxis(h(2),'left');
        h(2).YLim=[floor(min(data.Npas.tlim(Tint_6h).data)),ceil(max(data.Npas.tlim(Tint_6h).data))];
    end
    if ~isempty(data.B) && ~isempty(data.B.tlim(Tint_6h))
        yyaxis(h(2),'right');
        h(2).YLim=[floor(min(data.B.abs.tlim(Tint_6h).data)),ceil(max(data.B.abs.tlim(Tint_6h).data))];
    end
    irf_zoom(h(5:8),'y');
    
    %Remove overlapping ticks
    for iax=1:8
        cax=h(iax);
        mintick = min(cax.YTick);
        maxtick = max(cax.YTick);
        minlim = cax.YLim(1);
        maxlim = cax.YLim(2);
        
        if maxtick>0
            if maxlim<1.1*maxtick
                newmax = 1.1*maxtick;
            else
                newmax = maxlim;
            end
        else
            if abs(maxlim)>0.9*abs(maxtick)
                newmax = 0.9*maxtick;
            else
                newmax = maxlim;
            end
        end
        
        if mintick>0
            if minlim>0.9*mintick
                newmin = 0.9*mintick;
            else
                newmin = minlim;
            end
        else
            if abs(minlim)<1.1*abs(mintick)
                newmin=1.1*mintick;
            else
                newmin=minlim;
            end
        end
        cax.YLim=[newmin,newmax];
    end

    %Update text
    if ~isempty(data.solopos.tlim(Tint_6h))
        teststr = ['SolO: ',[sprintf('%.2f',data.solopos.tlim(Tint_6h).data(1,1)/Au),'Au, '],...
            [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,3)*180/pi)),'°, '],...
            [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,2)*180/pi)),'°']];        
        text1.String=teststr;
    else
        text1.String=[];
    end
    
    filesmth = Tint_6h(1);
    filesmth = filesmth.utc;
    filestr1 = filesmth(1:13);
    filestr1([5,8])=[];
    
    filesmth = Tint_6h(end);
    filesmth = filesmth.utc;
    filestr2 = filesmth(1:13);
    filestr2([5,8])=[];
    path2=fullfile(paths.path_6h,[filestr1,'_',filestr2,'.png']);
    print('-dpng',path2);
    
    %Print 2h figures
    for i2h=1:3
        %Define 2h interval and zoom in
        Tint_2h = Tint_6h(1)+[60*60*2*(i2h-1),60*60*2*(i2h)];
        irf_zoom(h(1:8),'x',Tint_2h);
        irf_zoom(h(1),'y');
        
        Neflag = ~isempty(data.Ne.tlim(Tint_2h));
        Npasflag = ~isempty(data.Npas) && ~isempty(data.Npas.tlim(Tint_2h));
        if Neflag && Npasflag
            yyaxis(h(2),'left');
            h(2).YLim=[min(floor([min(data.Npas.tlim(Tint_2h).data),min(data.Ne.tlim(Tint_2h).data)])),...
                max(ceil([max(data.Npas.tlim(Tint_2h).data),max(data.Ne.tlim(Tint_2h).data)]))];
        elseif Neflag
            yyaxis(h(2),'left');
            h(2).YLim=[floor(min(data.Ne.tlim(Tint_2h).data)),ceil(max(data.Ne.tlim(Tint_2h).data))];
        elseif Npasflag
            yyaxis(h(2),'left');
            h(2).YLim=[floor(min(data.Npas.tlim(Tint_2h).data)),ceil(max(data.Npas.tlim(Tint_2h).data))];
        end
        if ~isempty(data.B) && ~isempty(data.B.tlim(Tint_2h))
            yyaxis(h(2),'right');
            h(2).YLim=[floor(min(data.B.abs.tlim(Tint_2h).data)),ceil(max(data.B.abs.tlim(Tint_2h).data))];
        end
        
        irf_zoom(h(5:8),'y');
        %Remove overlapping Tics
        for iax=1:8
            cax=h(iax);
            mintick = min(cax.YTick);
            maxtick = max(cax.YTick);
            minlim = cax.YLim(1);
            maxlim = cax.YLim(2);
            
            if maxtick>0
                if maxlim<1.1*maxtick
                    newmax = 1.1*maxtick;
                else
                    newmax = maxlim;
                end
            else
                if abs(maxlim)>0.9*abs(maxtick)
                    newmax = 0.9*maxtick;
                else
                    newmax = maxlim;
                end
            end
            
            if mintick>0
                if minlim>0.9*mintick
                    newmin = 0.9*mintick;
                else
                    newmin = minlim;
                end
            else
                if abs(minlim)<1.1*abs(mintick)
                    newmin=1.1*mintick;
                else
                    newmin=minlim;
                end
            end
            cax.YLim=[newmin,newmax];
        end
        
        
        %Update text
        if ~isempty(data.solopos.tlim(Tint_2h))
            teststr = ['SolO: ',[sprintf('%.2f',data.solopos.tlim(Tint_6h).data(1,1)/Au),'Au, '],...
                [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,3)*180/pi)),'°, '],...
                [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,2)*180/pi)),'°']];
            text1.String=teststr;
        else
            text1.String=[];
        end
        
        filesmth = Tint_2h(1);
        filesmth = filesmth.utc;
        filestr1 = filesmth(1:13);
        filestr1([5,8])=[];
        
        filesmth = Tint_2h(end);
        filesmth = filesmth.utc;
        filestr2 = filesmth(1:13);
        filestr2([5,8])=[];
        path2=fullfile(paths.path_2h,[filestr1,'_',filestr2,'.png']);
        print('-dpng',path2);
    end
    
end
close(fig);


%%
