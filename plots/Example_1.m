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
tint=irf.tint('2006-09-27T17:17:00Z/2006-09-27T17:24:00Z'); % time interval

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1 % put to 0 if data already downloaded !!!!
    caa_download(tint,'C1_CP_FGM_5VPS')
    caa_download(tint,'C1_CP_CIS-HIA_ONBOARD_MOMENTS')
    caa_download(tint,'C1_CP_CIS_HIA_HS_1D_PEF')
    caa_download(tint,'C1_CP_RAP_ESPCT6')
    caa_download(tint,'C1_CP_PEA_PITCH_SPIN_DEFlux')
end

%%%%%%%%%%%%%%%%%%%%%%%
%% initialize figure
h=irf_plot(5,'newfigure'); % initialize new figure with 5 subplots

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('FGM B GSM');
% read data
B    = irf_get_data('B_vec_xyz_gse__C1_CP_FGM_5VPS','caa','ts');
gsmB = irf_gse2gsm(B);
% plot
irf_plot(hca,gsmB);
ylabel(hca,'B_{GSM} [nT] GSM');
irf_zoom(hca,'y',[-25 15])
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('CIS V');
V=irf_get_data('velocity_gse__C1_CP_CIS_HIA_ONBOARD_MOMENTS','caa','ts');
gsmV=irf_gse2gsm(V);
irf_plot(hca,gsmV)
ylabel(hca,'V [km/s] GSM');
irf_zoom(hca,'y',[-300 500])
irf_legend(hca,{'V_X','V_Y','V_Z'},[0.2 0.95])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('CIS spectrogram');
irf_plot(hca,'flux__C1_CP_CIS_HIA_HS_1D_PEF','colorbarlabel',{'log_{10} dEF','keV/cm^2 s sr keV'},'fitcolorbarlabel');
caxis([3.9 6.1]);
set(hca,'yscale','log')
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
ylabel(hca,'E [eV]')
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
irf_colormap('default');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('RAPID spectrogram');
irf_plot(hca,'Electron_Dif_flux__C1_CP_RAP_ESPCT6','colorbarlabel',{'log10 dF','1/cm^2 s sr keV'},'fitcolorbarlabel');
caxis([0.51 4.49]);
ylabel(hca,'E [keV]');
irf_legend(hca,{'C1'},[0.98 0.9],'color','k');
set(hca,'yscale','log');
set(hca,'ytick',[5e1 1e2 2e2 5e2 1e3])


%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('PEACE spectrogram');
irf_plot(hca,'Data__C1_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel',{'log10 dEF','keV/cm^2 s sr keV'},'fitcolorbarlabel');
caxis([5.9 7.6]);
set(hca,'yscale','log','ylim',[100 3e4])
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
ylabel(hca,'E [eV]');

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_plot_axis_align         % align the width of all panels
irf_zoom(h,'x',tint);       % zoom all panels to the same time interval 
irf_pl_number_subplots(h);  % add a), b), c) to panels
irf_timeaxis(h);            % add timeaxis ticksmarks and labels
irf_legend(h(1),'Example 1',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% add tmarks and mark intervals
% add line marks
tmarks=irf_time([2006 9 27 17 17 30;2006 9 27 17 18 20;2006 9 27 17 19 45;2006 9 27 17 21 0;2006 9 27 17 23 0]);
irf_pl_mark(h,tmarks,'black','LineWidth',0.5)
text_tmarks={'A','B','C','D','E'};
ypos=get(h(1),'ylim');ypos(1)=[];
for j=1:length(tmarks)
    irf_legend(h(1),text_tmarks{j},[tmarks(j) ypos],'horizontalalignment','center');
end
% add interval mark
tmarks=irf_time([2006 9 27 17 25 0])+[0 5*60];
irf_pl_mark(h(1:2),tmarks)


%%%%%%%%%%%%%%%%%%%%%%%%
% To print the figure 
%
% matlab> set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% to get bitmap file
% matlab> print -dpng delme.png
% to get pdf file with no white margins:
% 1) print eps file from matlab
% matlab> print -depsc2 -painters -loose delme.eps  % -loose is necessary for MATLAB2014 and later
% 2) from terminal convert to eps file without white margins
% epstool routine can be installed, e.g. on mac >brew install epstool
% terminal> epstool --copy --bbox delme.eps delme_crop.eps
% 3) convert eps file to pdf, result is in delme_crop.pdf
% terminal> ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None delme_crop.eps

%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary directory
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('When finnished with the example, ');
disp('remove the temporary directory in which you are located!')
disp('>p=pwd;cd ..; rmdir(p,''s'');');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')


