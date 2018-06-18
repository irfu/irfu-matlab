% you have to be located in a new directory and download data (see later)

% define time interval
tint=[toepoch([2006 9 27 17 27 0]) toepoch([2006 9 27 17 28 30])];

% to download data from CAA execute once the next line (uncommented)
% caa_download(tint,'*PEA*PITCH_*');

% loops "if 1,... end" are convenient when folding is enabled to fast browse the code 
if 1 % initialize figure
    close 61; fn=figure(61);
    set(fn,'Position',[10 400 500 700])
    h(1)=axes('position',[0.15 0.75 0.7 0.2]); % [x y dx dy]
    h(2)=axes('position',[0.15 0.51 0.7 0.2]); % [x y dx dy]
    h(3)=axes('position',[0.15 0.1 0.3 0.3]); % [x y dx dy]
    h(4)=axes('position',[0.62 0.1 0.3 0.3]); % [x y dx dy]
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',1);
end

if 1   % PANEL: PEACE PITCH_SPIN_DEFlux spectrogram omni
    hca=irf_panel('C3 PEACE energy spectra');
    varname='Data__C3_CP_PEA_PITCH_SPIN_PSD';
    [~,dobj,~,varunits]=c_caa_var_get(varname);
    varunits=['log10 PSD [' varunits ']'];
    plot(hca,dobj,varname,'sum_dim1','colorbarlabel',varunits,'fitcolorbarlabel');
    caxis(hca,[-7 2]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E [eV]');
end
if 1   % PANEL: PEACE PEA_PITCH_3DRH_PSD high res
    hca=irf_panel('C3 PEACE 3DR energy spectra');
    varname='Data__C3_CP_PEA_PITCH_3DRH_PSD';
    res=c_caa_construct_subspin_res_data(varname);
    [delmett,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['log10 PSD [' res.dataunits ']']);
    if 1 % energy spectorgram (integrated over polar angles)
        specrec.f=(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label=[res.enlabel];
        irf_spectrogram(hca,specrec);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E [eV]');
    elseif 0  %#ok<UNRCH> % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=15;
        res.en(enindex)
        specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        irf_spectrogram(hca,specrec);
        set(hca,'ytick',[30 60 90 120 150]);
    end
    caxis(hca,[-7 2]);
    irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
    set(hca,'yscale','log');
end
irf_colormap(hca,'default');
if 1 % plotting 1st distribution function 
    hca=irf_panel('fdistr');
    tt=irf_time([2006 9 27 17 27 03]) ; % time around which check for 3DR and other
    [p1,~,p1m]=c_caa_var_get('Data__C3_CP_PEA_PITCH_SPIN_PSD');
    xx=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_PSD');
    en=xx.data(1,:);ind_en1=find(en>0);en1=en(ind_en1);
    xx=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_3DRH_PSD');
    en=xx.data(1,:);ind_en2=find(en>0);en2=en(ind_en2);
    [p2,~,p2m]=c_caa_var_get('Data__C3_CP_PEA_PITCH_3DRH_PSD');
    ind1=find(abs(p1m.t-tt)<3); % find time within 3s of tt
    ind2=find(abs(p2m.t-tt)<3);
    for j=1:size(p2.data,2)
        for jj=1:size(p2.data,3)
            loglog(hca,en2,squeeze(p2.data(ind2,j,jj,ind_en2)),'r.');hold(hca,'on');
        end
    end
    for j=1:size(p1.data,2)
        loglog(hca,en1,squeeze(p1.data(ind1,j,ind_en1)),'b.');hold(hca,'on');
    end
    irf_legend(hca,'3DR',[0.9,0.98],'color','r');
    irf_legend(hca,'SPIN',[0.9,0.88],'color','b');
    grid(hca,'on');
    set(hca,'xlim',[30 30e3],'ylim',[1e-4 1e3]);
    set(hca,'xtick',[1 1e1 1e2 1e3 1e4 1e5]);
    xlabel(hca,'Energy [eV]');
    ylabel(hca,'PSD [s^3/km^6]');
end
if 1 % plotting 3DRH distribution functions
    hca=irf_panel('fdistr2');
    tt=irf_time([2006 9 27 17 28 20]) ; % time around which check for 3DR and other
    [p1,~,p1m]=c_caa_var_get('Data__C3_CP_PEA_PITCH_SPIN_PSD');
    xx=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_PSD');
    en=xx.data(1,:);ind_en1=find(en>0);en1=en(ind_en1);
    xx=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_3DRH_PSD');
    en=xx.data(1,:);ind_en2=find(en>0);en2=en(ind_en2);
    [p2,~,p2m]=c_caa_var_get('Data__C3_CP_PEA_PITCH_3DRH_PSD');
    ind1=find(abs(p1m.t-tt)<3);
    ind2=find(abs(p2m.t-tt)<3);
    for j=1:size(p2.data,2)
        for jj=1:size(p2.data,3)
            loglog(hca,en2,squeeze(p2.data(ind2,j,jj,ind_en2)),'r.');hold(hca,'on');
        end
    end
    for j=1:size(p1.data,2)
        loglog(hca,en1,squeeze(p1.data(ind1,j,ind_en1)),'b.');hold(hca,'on');
    end
    irf_legend(hca,'3DR',[0.9,0.98],'color','r');
    irf_legend(hca,'SPIN',[0.9,0.88],'color','b');
    grid(hca,'on');
    set(hca,'xlim',[30 30e3],'ylim',[1e-4 1e3]);
    set(hca,'xtick',[1 1e1 1e2 1e3 1e4 1e5]);
    xlabel(hca,'Energy [eV]');
    ylabel(hca,'PSD [s^3/km^6]');
end

% some general code to adjust plots, add labels ...
irf_plot_axis_align(1,h(1:2))
irf_zoom(h(1:2),'x',tint);
irf_pl_number_subplots(h,[0.02,1.01],'fontsize',14);
irf_timeaxis(h(1:2));
irf_legend(0,'Example 10',[0,0],'color',[0.5 0.5 0.5])
