
%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval
tint=[irf_time([2006 9 27 17 10 0]) irf_time([2006 9 27 17 40 0])]; % time interval

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once)
if 0, % put to 1 if needed to download
    caa_download(tint,'C1_CP_FGM_5VPS')
    caa_download(tint,'C1_CP_CIS-HIA_ONBOARD_MOMENTS')
    caa_download(tint,'C1_CP_CIS_HIA_HS_1D_PEF')
    caa_download(tint,'C1_CP_RAP_ESPCT6')
    caa_download(tint,'C1_CP_PEA_PITCH_SPIN_DEFlux')
    caa_download % repeat until all data are downloaded
    disp('Execute "caa_download" until all data are downloaded');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%
% initialize figure
h=irf_plot(5); % 5 subplots
i_subplot=1;

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=h(i_subplot); i_subplot=i_subplot+1;
% read data
[caaB,~,B]=c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_5VPS');
gsmB=irf_gse2gsm(B);
% plot
irf_plot(hca,gsmB);
ylabel(hca,'B [nT] GSM');
irf_zoom(hca,'y',[-25 15])
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.95])
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=h(i_subplot); i_subplot=i_subplot+1;
[caaV,~,V]=c_caa_var_get('velocity_gse__C1_CP_CIS_HIA_ONBOARD_MOMENTS');
gsmV=irf_gse2gsm(V);
irf_plot(hca,gsmV)
ylabel(hca,'V [km/s] GSM');
irf_zoom(hca,'y',[-300 500])
irf_legend(hca,{'V_X','V_Y','V_Z'},[0.02 0.1])
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=h(i_subplot); i_subplot=i_subplot+1;
irf_plot(hca,'flux__C1_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
caxis([3.9 6.1]);
set(hca,'yscale','log')
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
ylabel(hca,'E [eV]')
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
irf_colormap('default');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=h(i_subplot); i_subplot=i_subplot+1;
irf_plot(hca,'Electron_Dif_flux__C1_CP_RAP_ESPCT6','colorbarlabel','log10 dF\newline 1/cm^2 s sr keV','fitcolorbarlabel');
caxis([0.51 4.49]);
ylabel(hca,'E [keV]');
irf_legend(hca,{'C1'},[0.98 0.9],'color','k');
set(hca,'yscale','log');
set(hca,'ytick',[5e1 1e2 2e2 5e2 1e3])


%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=h(i_subplot); i_subplot=i_subplot+1;
irf_plot(hca,'Data__C1_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
caxis([5.9 7.6]);
set(hca,'yscale','log','ylim',[100 3e4])
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
ylabel('E [eV]');

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_pl_number_subplots(h);
irf_timeaxis(h);
irf_legend(h(1),'Example 1',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% add tmarks and mark intervals
% add line marks
tmarks=irf_time([2006 9 27 17 12 0;2006 9 27 17 15 0;2006 9 27 17 18 0;2006 9 27 17 21 0;2006 9 27 17 23 0]);
irf_pl_mark(h,tmarks,'black','LineWidth',0.5)
text_tmarks={'A','B','C','D','E'};
ypos=ylim(h(1));ypos(2)=ypos(2);ypos(1)=[];
for j=1:length(tmarks)
    irf_legend(h(1),text_tmarks{j},[tmarks(j) ypos],'horizontalalignment','center');
end
% add interval mark
tmarks=irf_time([2006 9 27 17 25 0])+[0 5*60];
irf_pl_mark(h(1:2),tmarks)


%%%%%%%%%%%%%%%%%%%%%%%%
% print the figure
if 0, % put to 1 if you want to print
    set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
    print -dpng -painters Example_1.png;
end