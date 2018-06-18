%%%%%%%%%%%%%%%%%%%%%%%
% go to new/empty directory 
% >cd new_directory
% here using temporary directory
tempdir_name=tempname;
mkdir(tempdir_name);
cd(tempdir_name);
disp(['Moving to temporary directory: ' tempdir_name]); 

%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval
tint=irf_time([2001 2 26 5 18 0]) + [0 60];

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1 % put to 0 if data already downloaded !!!!
    caa_download(tint,'C1_CP_FGM_FULL_ISR2')
    caa_download(tint,'C1_CP_WHI_NATURAL')
    caa_download(tint,'C1_CP_EFW_L2_E','nowildcard')
    caa_download(tint,'C1_CP_STA_CWF_HBR_ISR2')
    caa_download(tint,'C1_CP_CIS_HIA_HS_1D_PEF')
    caa_download(tint,'C1_CP_PEA_PITCH_SPIN_DEFlux')
    download_status=caa_download; % repeat until all data are downloaded
    if download_status==0 % some data are still in queue
      disp('___________!!!!_____________')
      disp('Some data where put in queue!')
      disp('To see when they are ready and to download execute "caa_download".');
      return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% initialize figure
h=irf_plot(7,'newfigure'); % 5 subplots

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('FGM B ISR2');
% read data
B=irf_get_data('B_vec_xyz_isr2__C1_CP_FGM_FULL_ISR2','caa','mat');
% plot
irf_plot(hca,B);
ylabel(hca,{'B FGM','[nT] ISR2'});
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('EFW E ISR2 c?');
diE  = irf_get_data('E_Vec_xy_ISR2__C1_CP_EFW_L2_E' ,'caa','mat');
irf_plot(hca,diE(:,1:3));
ylabel(hca,{'E EFW','[mV/m] ISR2'});
irf_zoom(hca,'y');
irf_legend(hca,{'E_X','E_Y'},[0.98 0.05])
irf_legend(hca,{'C1'},[0.02 0.95],'color','k')

% %%%%%%%%%%%%%%%%%%%%%%%
% % new panel
 hca=irf_panel('STAFF B x');
% % read data
Bstaff=irf_get_data('B_vec_xyz_Instrument__C1_CP_STA_CWF_HBR_ISR2','caa','mat');
fmin=1;
Bstaff=irf_filt(Bstaff,fmin,0,[],5);
res=irf_ebsp([],Bstaff,[],B,[],[1 180],'polarization','fac');
specrec = struct();
specrec.t = res.t;
specrec.f = res.f;
specrec.p = res.ellipticity;specrec.p_label='ellipticity';
specrec.p(res.dop < 0.7) = NaN;
%specrec.p = res.dop;specrec.p_label=' ';
%specrec.p = res.dop2d;specrec.p_label=' ';
%specrec.p = res.planarity;specrec.p_label=' ';
%specrec.p = squeeze(res.bb(:,:,1));specrec.p_label = 'nT^2/Hz';
% % plot
irf_spectrogram(hca,specrec,'lin');
irf_colormap(hca,'poynting');
irf_zoom(hca,'y',[0 200]);
% ylabel(hca,{'B STAFF','[nT] ISR2'});
% irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
% irf_legend(hca,{'C1'},[0.98 0.98],'color','k')
% irf_legend(hca,{'high pass filter at' num2str(fmin,'%3.1f') ' Hz'},[0.02 0.02],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('STAFF B spectra');
% read data
Bstaff=irf_get_data('B_vec_xyz_Instrument__C1_CP_STA_CWF_HBR_ISR2','caa','mat');
fmin=1;
Bstaff=irf_filt(Bstaff,fmin,0,[],5);
% construct spectra
specrec = irf_powerfft(Bstaff(:,[1 4]),256,450);
specrec.p_label ='[nT]^2/Hz';
% plot
irf_spectrogram(hca,specrec);
irf_zoom(hca,'y',[0 200]);
caxis(hca,[-7.9 -4.1]);
ylabel(hca,{'B STAFF','[nT] ISR2'});
irf_legend(hca,{'B_Z'},[0.98 0.1])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')
irf_colormap(hca,'default');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('WHISPER spectrogram natural');
varname='Electric_Spectral_Power_Density__C1_CP_WHI_NATURAL';
varunits='(V/m)^2/Hz';
% REMOVE 'fillspectgrogramgaps' flag in the next line if precise intervals of
% WHISPER measurements are needed !!!!
% If working with shorter intervals can also remove 'tint' option
irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','fillspectrogramgaps','nolabels');
ylabel(hca,{'E WHISPER','[mV/m] ISR2'});
% polish the panel
caxis(hca,[-15.9 -11.1]);
set(hca,'yscale','log','ytick',[3 4 5 1e1 20 30 50 ]);
irf_zoom(hca,'y',[2 12]);


%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('CIS spectrogram');
irf_plot(hca,'flux__C1_CP_CIS_HIA_HS_1D_PEF','colorbarlabel',{'log_{10} dEF','keV/cm^2 s sr keV'},'fitcolorbarlabel');
caxis([3.9 6.1]);
set(hca,'yscale','log')
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
ylabel(hca,{'E ions CIS/HIA','[eV]'})
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
irf_colormap('default');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('PEACE spectrogram');
irf_plot(hca,'Data__C1_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel',{'log10 dEF','keV/cm^2 s sr keV'},'fitcolorbarlabel');
caxis([5.9 7.6]);
set(hca,'yscale','log','ylim',[100 3e4])
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
ylabel(hca,{'E e- PEACE','[eV]'});

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_pl_number_subplots(h);
irf_timeaxis(h);
irf_legend(h(1),'Example 4',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% to print the figure uncomment the lines below
%
% set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% print -dpng -painters Example_4.png;

%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary directory
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('When finnished with the example, ');
disp('remove the temporary directory in which you are located!')
disp('>p=pwd;cd ..; rmdir(p,''s'');');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')


