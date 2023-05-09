function quicklooks_24_6_2_h(data,paths,Tint,logoPath)
% Given data in the struct 'data' (see solo.qli.quicklooks_main), generates
% plots and saves them in the paths specified in the struct 'paths' (see
% solo.qli.quicklooks_main). Computes spectrum of B, so takes a while to run.
% Tint should be a 24hour time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-02T00:00:00.00Z');

tBeginSec = tic();



% Setup figure:
lwidth=1.0;
fsize=18;
legsize=22;
h=irf_plot(10,'newfigure');
fig=gcf;
fig.Position=[1,1,1095,800];
colors = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];

Units=irf_units;

Me=Units.me; %Electron mass [kg]
epso = Units.eps0; %Permitivitty of free space [Fm^-1]
mp=Units.mp; %proton mass [km]
qe=Units.e; %elementary charge [C]

%==============
% Fill panel 1
%==============
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint),'linewidth',lwidth);
    hold(h(1),'on');
    irf_plot(h(1),data.B.abs.tlim(Tint),'linewidth',lwidth);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',legsize);
ylabel(h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',fsize);

tBeginSec = log_time('End panel 1', tBeginSec);



%==============
% Fill panel 2
%==============
%%
hold(h(2),'on');
if ~isempty(data.Ne)
    irf_plot(h(2),data.Ne.tlim(Tint),'-','color',colors(1,:),'linewidth',lwidth);
end
if ~isempty(data.Npas)
    irf_plot(h(2),data.Npas.tlim(Tint),'-','color',colors(2,:),'linewidth',lwidth);
end
ylabel(h(2),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',fsize);
h(2).ColorOrder=colors;
irf_legend(h(2),{'N_{e,RPW}','N_{i,PAS}','|B|'},[0.98 0.16],'Fontsize',legsize);

yyaxis(h(2),'right');
if ~isempty(data.B)
    fci = qe*data.B.abs*10^-9/mp/(2*pi);
    irf_plot(h(2),data.B.abs.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
    %Bnan = rmmissing(data.B.abs.data);
    %if ~isempty(Bnan)
    %    h(2).YLim=[floor(min(abs(Bnan))),ceil(max(abs(Bnan)))];
    %end
    minAbsB = min(data.B.tlim(Tint).abs.data);
    maxAbsB = max(data.B.tlim(Tint).abs.data);
    if ~isnan(minAbsB) && ~isnan(maxAbsB)
        % Only zoom if min & max are not NaN (==> Avoid crash).
        irf_zoom(h(2),'y',[minAbsB-1, maxAbsB+1]);
    end
end
ylabel(h(2),{'|B|';'(nT)'},'interpreter','tex','fontsize',fsize);
h(2).YColor=[1,0,0];

tBeginSec = log_time('End panel 2', tBeginSec);



%======================================
% Fill panel 3 & 4: Spectra h(3), h(4)
%======================================
%%
if ~isempty(data.B)
   if  ~isempty(rmmissing(data.B.data))
    bb = data.B;
    if median(diff((bb.time.epochUnix))) < 0.1250*0.95
        fMag = 128; fMax = 7;
    else
        fMag = 8; fMax = 3;
    end
    b0 = bb.filt(0, 0.01,fMag, 5);

    % IMPORTANT NOTE: The call to irf_ebsp() seems very slow.
    tBeginSec = log_time('irf_ebsp(): Begin call', tBeginSec);
    ebsp = irf_ebsp([],bb,[],b0,[],[0.05 fMax],'fullB=dB', 'polarization', 'fac');
    tBeginSec = log_time('irf_ebsp(): End call', tBeginSec);

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
end

tBeginSec = log_time('End panel 3 & 4', tBeginSec);



%===============================
% Fill panel 5: Ion temperature
%===============================
if ~isempty(data.Tpas)
    irf_plot(h(5),data.Tpas.tlim(Tint),'color',colors(2,:),'linewidth',lwidth); 
end
irf_zoom(h(5),'y');
ylabel(h(5),{'T_i';'(eV)'},'interpreter','tex','fontsize',fsize);

tBeginSec = log_time('End panel 5', tBeginSec);



%==================================
% Fill panel 6: y,z PAS velocities
%==================================
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.y.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
    hold(h(6),'on');
    irf_plot(h(6),data.Vpas.z.tlim(Tint),'color',colors(3,:),'linewidth',lwidth);
end
irf_legend(h(6),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(6),'y');
ylabel(h(6),{'V_{T,N}';'(km/s)'},'interpreter','tex','fontsize',fsize);

tBeginSec = log_time('End panel 6', tBeginSec);



%=====================================
% Fill panel 7: Vrpw, Vpas velocities
%=====================================
hold(h(7),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(7),-data.Vrpw,'o-','color',colors(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(7),data.Vpas.x.tlim(Tint),'color',colors(2,:),'linewidth',lwidth);
end
irf_legend(h(7),{'V_{RPW}','V_{PAS}'},[0.98 0.18],'Fontsize',legsize);
irf_zoom(h(7),'y');
ylabel(h(7),{'V_R';'(km/s)'},'interpreter','tex','fontsize',fsize);

tBeginSec = log_time('End panel 7', tBeginSec);



%==============================
% Fill panel 8: Electric field
%==============================
if ~isempty(data.E)
    irf_plot(h(8),data.E.y,'color',colors(2,:),'linewidth',lwidth)
    hold(h(8),'on');
    %irf_plot(h(8),data.E.z,'color',colors(3,:),'linewidth',lwidth)
end
irf_legend(h(8),{'','E_y'},[0.98 0.20],'Fontsize',legsize);
irf_zoom(h(8),'y');
ylabel(h(8),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',fsize);

tBeginSec = log_time('End panel 8', tBeginSec);



%===================================
% Fill panel 9: Ion energy spectrum
%===================================
if ~isempty(data.ieflux)
    myFile=solo.db_list_files('solo_L2_swa-pas-eflux',Tint);
    iDEF   = struct('t',  data.ieflux.tlim(Tint).time.epochUnix);
    % for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
    for ii = 1:length(myFile)
        iEnergy = cdfread([myFile(ii).path '/' myFile(ii).name],'variables','Energy');
        iEnergy = iEnergy{1};
        iDEF.p = data.ieflux.data;
          
    end
    iDEF.p_label={'dEF','keV/','(cm^2 s sr keV)'};
    iDEF.f = repmat(iEnergy,1,numel(iDEF.t))';
    irf_spectrogram(h(9),iDEF,'log','donotfitcolorbarlabel');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    %caxis(h(9),[-1 1])
    hold(h(9),'on');
    set(h(9), 'YScale', 'log');
    colormap(h(9),jet)
    ylabel(h(9),{'W_{i}';'(eV)'},'interpreter','tex','fontsize',fsize);
end

tBeginSec = log_time('End panel 9', tBeginSec);



%=======================================
% Fill panel 10: E-field spectrum (TNR)
%=======================================
if ~isempty(data.Etnr)
    try
        [TNR] = solo.read_TNR(Tint);
    catch Exc
        if strcmp(Exc.identifier, 'read_TNR:FileNotFound')
            TNR = [];
        end
    end
    if isa(TNR,'struct')
        sz_tnr = size(TNR.p);
        if sz_tnr(1) == length(TNR.t) && sz_tnr(2) == length(TNR.f)
            irf_spectrogram(h(10),TNR,'log','donotfitcolorbarlabel')
            hold(h(10),'on');
            if ~isempty(data.Ne)
                %Electron plasma frequency
                wpe_sc = (sqrt(((data.Ne.tlim(Tint)*1000000)*qe^2)/(Me*epso)));                         
                fpe_sc = (wpe_sc/2/pi)/1000;
                irf_plot(h(10),fpe_sc,'r','linewidth',lwidth);
                fpe_sc.units = 'kHz';
                fpe_sc.name = 'f [kHz]';
            end
            hold(h(10),'off');
            text(h(10),0.01,0.3,'f_{pe,RPW}','units','normalized','fontsize',18,'Color','r');
            set(h(10), 'YScale', 'log'); 
            %set(h(10),'ColorScale','log')
            %caxis(h(10),[.01 1]*10^-12)
            ylabel(h(10),{'f';'(kHz)'},'interpreter','tex','fontsize',fsize);
            colormap(h(10),jet)  
            yticks(h(10),[10^1 10^2]);
            irf_zoom(h(10),'y',[10^1 10^2])
        end
    end
end

if isempty(data.Vrpw) && isempty(data.E) && isempty(data.Ne) && isempty(data.B) ...
        && isempty(data.Tpas) && isempty(data.Npas) && isempty(data.ieflux) ...
        && isempty(data.Etnr)
    nanPlot = irf.ts_scalar(Tint,ones(1,2)*NaN);
    irf_plot(h(10),nanPlot);
    grid(h(10),'off');
    ylabel(h(10),{'f';'(kHz)'},'interpreter','tex','fontsize',fsize);
end 

tBeginSec = log_time('End panel 10', tBeginSec);



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:10));
irf_zoom(h(1:10),'x',Tint);
irf_zoom(h(1),'y');
irf_zoom(h(5:10),'y');

h(2).YLabel.Position=[1.05,0.5,0];
yyaxis(h(2),'left');
h(2).YLabel.Units='normalized';
h(2).YLabel.Position=h(3).YLabel.Position;
% Add spacecraft position as text.
Au=149597871; %km
if ~isempty(data.solopos.tlim(Tint))
    teststr = ['SolO: ',[sprintf('%.2f',data.solopos.tlim(Tint).data(1,1)/Au),'Au, '],...
        [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint).data(1,3)*180/pi)),'\circ, '],...
        [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint).data(1,2)*180/pi)),'\circ']];
    text1=text(h(10),-0.11,-0.575,teststr,'units','normalized','fontsize',18);
else
    teststr=char();
    text1=text(h(10),-0.11,-0.575,teststr,'units','normalized','fontsize',18);
end

% Add Earth longitude as text.
if ~isempty(data.earthpos)
    teststr =['Earth: EcLon ',sprintf('%d',round(data.earthpos(1,2)*180/pi)),'\circ'];
    text2=text(h(10),-0.11,-0.925,teststr,'units','normalized','fontsize',18);
else
    teststr=char();
    text2=text(h(10),0.9,-0.925,teststr,'units','normalized','fontsize',18);
end

xtickangle(h(10),0)
% Add plot information and IRF logo
logopos = h(1).Position;
logopos(1)=logopos(1)+logopos(3)+0.06;
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

% Fix YTicks
for iax=1:10
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
%h(2).YLim=[0.8,200];

yyaxis(h(2),'right');
oldlims2_r=h(2).YLim;
oldticks2_r = h(2).YTick;
h(2).YScale='log';
h(2).YTick=[1,10,100];
%h(2).YLim=[0.1,200];

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
%h(2).YLim=oldlims2_r;
yyaxis(h(2),'left');
h(2).YScale='lin';
%h(2).YLim=oldlims2;
h(2).YTick=oldticks2;

h(5).YScale='lin';
h(5).YLim=oldlims5;
h(5).YTick=oldticks5;



%===========================
% Iterate over 6h intervals
%===========================
% Print 6h figures.
tBeginSec = log_time('Begin iterating over 6 h intervals', tBeginSec);
for i6h = 1:4
    
    %Zoom in to 6h interval and save plot.
    Tint_6h = Tint(1)+[60*60*6*(i6h-1),60*60*6*(i6h)];
    irf_zoom(h(1:10),'x',Tint_6h);
    irf_zoom(h(1),'y');
    %Zoom on N/|B| plot..
    Neflag   = ~isempty(data.Ne)   && ~isempty(data.Ne.tlim(Tint_6h)) && ~all(isnan(data.Ne.tlim(Tint_6h).data));
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
    if ~isempty(data.B) && ~isempty(data.B.tlim(Tint_6h)) && ~all(isnan(data.B.abs.tlim(Tint_6h).data))
        yyaxis(h(2),'right');
        h(2).YLim=[floor(min(data.B.abs.tlim(Tint_6h).data)),ceil(max(data.B.abs.tlim(Tint_6h).data))];
    end
    if ~isempty(data.Tpas.tlim(Tint_6h))
        minTi = min(data.Tpas.tlim(Tint_6h).abs.data);
        maxTi = max(data.Tpas.tlim(Tint_6h).abs.data);
        if ~isnan(minTi) && ~isnan(maxTi)
            % Only zoom if min & max are not NaN (==> Avoid crash).
            irf_zoom(h(5),'y',[minTi-2, maxTi+2]);
        end
    end
    if ~isempty(data.Vpas.tlim(Tint_6h))
        minVy = min(rmmissing(data.Vpas.y.tlim(Tint_6h).data));
        minVz = min(rmmissing(data.Vpas.z.tlim(Tint_6h).data));
        maxVy = max(rmmissing(data.Vpas.y.tlim(Tint_6h).data));
        maxVz = max(rmmissing(data.Vpas.z.tlim(Tint_6h).data));
        maxV = max(maxVy,maxVz);
        minV = min(minVy,minVz);
        if ~isempty(minV) && ~isempty(maxV)
            % Only zoom if min & max are not NaN (==> Avoid crash).
            irf_zoom(h(6),'y',[minV-10, maxV+10]);
        end
    end
    irf_zoom(h(7:8),'y');
    
    %Remove overlapping ticks
    for iax=1:10
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
            [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,3)*180/pi)),'\circ, '],...
            [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,2)*180/pi)),'\circ']];
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
    
    %==================================================
    % Iterate over 2h intervals within one 6h interval
    %==================================================
    tBeginSec = log_time('Begin iterating over 2 h intervals', tBeginSec);
    % Print 2h figures
    for i2h=1:3
        %Define 2h interval and zoom in
        Tint_2h = Tint_6h(1)+[60*60*2*(i2h-1),60*60*2*(i2h)];
        irf_zoom(h(1:10),'x',Tint_2h);
        irf_zoom(h(1),'y');
        
        Neflag   = ~isempty(data.Ne)   && ~isempty(data.Ne.tlim(Tint_2h)) && ~all(isnan(data.Ne.tlim(Tint_2h).data));
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
        if ~isempty(data.B) && ~isempty(data.B.tlim(Tint_2h)) && ~all(isnan(data.B.abs.tlim(Tint_2h).data))
            yyaxis(h(2),'right');
            h(2).YLim=[floor(min(data.B.abs.tlim(Tint_2h).data)),ceil(max(data.B.abs.tlim(Tint_2h).data))];
        end
        if ~isempty(data.Tpas.tlim(Tint_2h))
            minTi = min(data.Tpas.tlim(Tint_2h).abs.data);
            maxTi = max(data.Tpas.tlim(Tint_2h).abs.data);
            if ~isnan(minTi) && ~isnan(maxTi)
                % Only zoom if min & max are not NaN (==> Avoid crash).
                irf_zoom(h(5),'y',[minTi-2, maxTi+2]);
            end
        end
        if ~isempty(data.Vpas.tlim(Tint_2h))
            minVy = min(rmmissing(data.Vpas.y.tlim(Tint_2h).data));
            minVz = min(rmmissing(data.Vpas.z.tlim(Tint_2h).data));
            maxVy = max(rmmissing(data.Vpas.y.tlim(Tint_2h).data));
            maxVz = max(rmmissing(data.Vpas.z.tlim(Tint_2h).data));
            maxV = max(maxVy,maxVz);
            minV = min(minVy,minVz);
            if ~isempty(minV) && ~isempty(maxV)
                % Only zoom if min & max are not NaN (==> Avoid crash).
                irf_zoom(h(6),'y',[minV-10, maxV+10]);
            end
        end
        irf_zoom(h(7:8),'y');
        
        
        %Remove overlapping Tics
        for iax=1:10
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
                [' EcLat ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,3)*180/pi)),'\circ, '],...
                [' EcLon ',sprintf('%d',round(data.solopos.tlim(Tint_6h).data(1,2)*180/pi)),'\circ']];
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



tBeginSec = log_time('End of quicklooks_24_6_2_h.m', tBeginSec);


end



function tBeginSec = log_time(locationStr, tBeginSec)
    % Simple function for logging number of seconds from previous call.
    % For debugging speed.
    tSec = toc(tBeginSec);
    fprintf(1, '%s: %.1f [s]\n', locationStr, tSec)
    tBeginSec = tic();
end
