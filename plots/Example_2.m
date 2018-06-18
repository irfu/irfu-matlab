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
tint=[irf_time([2006 9 27 17 22 0]) irf_time([2006 9 27 17 24 30])];

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1 % put to 0 if data already downloaded !!!!
    caa_download(tint,'C?_CP_STA_PSD')
    caa_download(tint,'C?_CP_FGM_5VPS')
    caa_download(tint,'C?_CP_RAP_PAD_L3DD')
    download_status=caa_download; % repeat until all data are downloaded
    if download_status==0 % some data are still in queue
      disp('___________!!!!_____________')
      disp('Some data where put in queue!')
      disp('To see when they are ready and to download execute "caa_download".');
      return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%% initialize figure
h=irf_plot(6); % 5 subplots
irf_colormap('space');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('C1..C4 FGM BZ GSM');
% read data
c_eval('B?=irf_get_data(tint,''B_vec_xyz_gse__C?_CP_FGM_5VPS'',''caa'',''mat'');');
c_eval('gsmB?=irf_gse2gsm(B?);');
% plot
c_pl_tx(hca,'gsmB?',4)
ylabel(hca,{'B_Z GSM','[nT]'});
irf_zoom(hca,'y',[-1 17]);
irf_legend(hca,{'C1','C2','C3','C4'},[0.98, 0.95],'color','cluster');

%%%%%%%%%%%%%%%%%%%%%%%
% 4 panels of RAPID
for ic=1:4
  hca=irf_panel(['RAPID anisotropy C' num2str(ic)]);
  varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
  varunits='[1/cm^2 s sr keV]';
  irf_plot(varname,'ax',hca,'colorbarlabel',{'log_{10}dPF',varunits},'fitcolorbarlabel','comp_dim1','comp',2,'nolabels');
  caxis(hca,[1.1 3.9]);
  irf_legend(hca,['C' num2str(ic)],[0.98 0.98]);
  set(hca,'ytick',[45 90 135]);
  ylabel(hca,{'pitch angle','[deg]'});
end

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('STAFF B spectrogram');
ic=4;
varname=irf_ssub('BB_xxyyzz_sr2__C?_CP_STA_PSD',ic);
varunits='nT^2/Hz';
irf_plot(varname,'ax',hca,'colorbarlabel',{'log_{10}S',varunits},'fitcolorbarlabel','comp',1,'nolabels');
ylabel(hca,{'frequency','[Hz]'});
hold(hca,'on');
c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic);
irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
caxis(hca,[-9 -5]);
set(hca,'yscale','log','ylim',[11 599]);

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_pl_number_subplots(h);
irf_timeaxis(h);
irf_legend(0,'Example 2',[0.02 0.02],'color',[0.5 0.5 0.5]);


%%%%%%%%%%%%%%%%%%%%%%%%
% to print the figure uncomment the lines below
%
% set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% print -dpng -painters Example_2.png;

%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary directory
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('When finnished with the example, ');
disp('remove the temporary directory in which you are located!')
disp('>p=pwd;cd ..; rmdir(p,''s'');');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')

