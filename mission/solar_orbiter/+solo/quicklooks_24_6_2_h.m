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
h=irf_plot(9,'newfigure');
fig=gcf;
fig.Position=[1,1,1095,800];
colors = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),'B_{SRF} (nT)','interpreter','tex','fontsize',fsize);

if ~isempty(data.B)
    irf_plot(h(2),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
ylabel(h(2),'|B| (nT)','interpreter','tex','fontsize',fsize);

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
    ylabel(h(3),'f (Hz)','fontsize',fsize);
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
    ylabel(h(4),'f (Hz)','fontsize',fsize);
    
    crr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
    cgg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
    cbb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
    bgrcmap = [crr' cgg' cbb'];
    colormap(h(4),bgrcmap);
end


% Densities
hold(h(5),'on');
if ~isempty(data.Ne)
    irf_plot(h(5),data.Ne.tlim(Tint),'color',colors(1,:),'linewidth',lwidth);
else
end
if ~isempty(data.Npas)
    irf_plot(h(5),data.Npas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(5),'N (cm^{-3})','interpreter','tex','fontsize',fsize);
irf_legend(h(5),{'N_{e,RPW} ',' N_{i,PAS}'},[0.98 0.16],'Fontsize',legsize);
irf_zoom(h(5),'y');

if ~isempty(data.Tpas)
    irf_plot(h(6),data.Tpas.tlim(Tint),'color',colors(2,:),'linewdith',lwidth);
end
ylabel(h(6),'T_i (eV)','interpreter','tex','fontsize',fsize);
irf_zoom(h(6),'y');

% y,z PAS velocities
if ~isempty(data.Vpas)
    irf_plot(h(7),V_pas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(7),'on');
    irf_plot(h(7),V_pas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(7),{'','v_{y}','v_{z}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),'v_{yz} (km/s)','interpreter','tex','fontsize',fsize);

hold(h(8),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(8),data.Vrpw,'o-','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(8),V_pas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(8),{'V_{RPW}','V_{PAS}'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(8),'y');
ylabel(h(8),'v_{x} (km/s)','interpreter','tex','fontsize',fsize);

if ~isempty(data.E)
    irf_plot(h(9),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(9),'on');
    irf_plot(h(9),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(9),{'','E_y','E_z'},[0.98 0.15],'Fontsize',legsize);
irf_zoom(h(9),'y');
ylabel(h(9),'E (mV/m)','interpreter','tex','fontsize',fsize);



irf_plot_axis_align(h(1:9));
irf_zoom(h(1:9),'x',Tint);
irf_zoom(h(1:2),'y');
irf_zoom(h(5:9),'y');

% Add spacecraft position as text.
Au=149597871; %km
text1=text(h(9),-0.11,-0.5,['R=',sprintf('%.2f',data.solopos.tlim(Tint).data(1,1)/Au),'Au'],'units','normalized','fontsize',18);
text2=text(h(9),0.98,-0.5,['R=',sprintf('%.2f',data.solopos.tlim(Tint).data(end,1)/Au),'Au'],'units','normalized','fontsize',18);

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

% Print 6h figures
for i6h = 1:4
    
    %Zoom in to 6h interval and save plot.
    Tint_6h = Tint(1)+[60*60*6*(i6h-1),60*60*6*(i6h)];
    irf_zoom(h(1:9),'x',Tint_6h);
    irf_zoom(h(1:2),'y');
    irf_zoom(h(5:9),'y');
    
    %Update text
    if ~isempty(data.solopos.tlim(Tint_6h))
        text1.String=['R=',sprintf('%.2f',data.solopos.tlim(Tint_6h).data(1,1)/Au),'Au'];
        text2.String=['R=',sprintf('%.2f',data.solopos.tlim(Tint_6h).data(end,1)/Au),'Au'];
    else
        text1.String=[];
        text2.String=[];
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
        irf_zoom(h(1:9),'x',Tint_2h);
        irf_zoom(h(1:2),'y');
        irf_zoom(h(5:9),'y');
        
        %Update text
        if ~isempty(data.solopos.tlim(Tint_2h))
            text1.String=['R=',sprintf('%.2f',data.solopos.tlim(Tint_2h).data(1,1)/Au),'Au'];
            text2.String=['R=',sprintf('%.2f',data.solopos.tlim(Tint_2h).data(end,1)/Au),'Au'];
        else
            text1.String=[];
            text2.String=[];
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
