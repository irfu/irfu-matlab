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
tintIso='2002-08-21T07:50:00Z/2002-08-21T08:00:00Z';
tint = irf_time(tintIso,'utc>tint');
% time interval can be specified directly, uncomment below
%tint=[irf_time([2002 8 21 7 50 0]) irf_time([2002 8 21 8 0 0])]; 

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA or CSA (needed only once!!!!!)
if 1 % put to 0 if data already downloaded !!!!
	dataArchive = 'csa'; % can be 'caa' 
    caa_download(tint,'C1_CP_FGM_5VPS'             ,dataArchive)
    caa_download(tint,'C1_CP_PEA_PITCH_SPIN_DEFlux',dataArchive)
	caa_download(tint,'C1_CP_EFW_L3_P'             ,dataArchive)
    caa_download(tint,'C1_CP_STA_PSD'              ,dataArchive)
    caa_download(tint,'C1_CP_AUX_POSGSE_1M'        ,dataArchive)
end

%%%%%%%%%%%%%%%%%%%%%%%
% initialize figure
h=irf_plot(3); % 3 subplots

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('FGM B GSM');
% read data
B=irf_get_data('B_vec_xyz_gse__C1_CP_FGM_5VPS','caa','mat');
gsmB=irf_gse2gsm(B);
% plot
irf_plot(hca,gsmB(:,[1 4])); % select only Bz component
hold(hca,'on');
irf_plot(hca,[gsmB(1,1) 0;gsmB(end,1) 0],'b:'); % draw line at Y=0
set(hca,'xgrid','off','ygrid','off');
ylabel(hca,'B_Z [nT] GSM');
irf_zoom(hca,'y',[-25 15])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('PEACE spectrogram');
% read data
Vps = irf_get_data('Spacecraft_potential__C1_CP_EFW_L3_P','caa','mat');
Vsat = [Vps(:,1) -Vps(:,2)]; % satellite potential is negative of 'probe ti sc' potential 
% plot
irf_plot(hca,'Data__C1_CP_PEA_PITCH_SPIN_DEFlux',...
	'sum_dim1pitch',...% make average over 1 dimension, pitch angle average is not the same as simple average
	'colorbarlabel',{'log10 dEF','keV/cm^2 s sr keV'});
caxis(hca,[5.1 7.9]);
set(hca,'yscale','log','ylim',[10 3e4])
set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
hold(hca,'on');
irf_plot(hca,Vsat,'w-');
irf_legend(hca,{'C1'},[0.98 0.05],'color','k')
ylabel(hca,'E [eV]');

%%%%%%%%%%%%%%%%%%%%%%%
% new panel
hca=irf_panel('STAFF spectrogram B and fce/flh lines');
irf_plot(hca,'BB_xxyyzz_sr2__C1_CP_STA_PSD',...
	'colorbarlabel',{'log10 S',irf_get_data(varname,'caa','units')},...
	'sum_dim2'); % sum Bx^2, By^2, Bz^2
hold(hca,'on');
fce=irf_plasma_calc(B,0,0,0,0,'Fce'); % calculate electron gyrofrequency
irf_plot(hca,fce,'-','linewidth',0.2,'color','w');
flh=irf_plasma_calc(B,1,0,0,0,'Flh'); % calculate lower hybrid frequency (underdense case in outer magnetosphere)
irf_legend(hca,'f_{ce}',[0.02 0.9],'color','w');
irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
irf_legend(hca,'f_{lh}',[0.02 0.1],'color','k');
hold(hca,'off');
caxis(hca,[-10 -4]);
set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
irf_zoom(hca,'y',[10 4000]);

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_pl_number_subplots(h);
irf_timeaxis(h);
irf_legend(h(1),'Example CSA',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% add position markers to x labels, remove UT date at bottom and add at top
R=irf_get_data('sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','caa','mat');
irf_timeaxis(h(end),'usefig',[R(:,1) R(:,2:4)/6371.2],...
	{'X [Re] GSE','Y [Re] GSE','Z [Re] GSE'})
irf_timeaxis(h(end),'nodate')
title(h(1),irf_time(tint(1),'yyyy-mm-dd'));

%%%%%%%%%%%%%%%%%%%%%%%%
% add tmarks and mark intervals
% add line marks
tmarks=irf_time([2002 8 21 7 53 15;2002 8 21 7 53 55]);
irf_pl_mark(h,tmarks,'black','LineWidth',0.5)
text_tmarks={'A','B','C','D','E'};
ypos=ylim(h(1));ypos(2)=ypos(2);ypos(1)=[];
for j=1:length(tmarks)
    irf_legend(h(1),text_tmarks{j},[tmarks(j) ypos],'horizontalalignment','center');
end


%%%%%%%%%%%%%%%%%%%%%%%%
% to print the figure uncomment appropriate lines below
%
% set(gcf,'paperpositionmode','auto')  % to get the same on paper as on screen
% print -dpng -painters Example_1.png; % to print png file
% print -depsc2 -painters delme.eps; % to print eps file with no white margins
%	to obtain pdf without white margins execute on the system command one of:
% ps2pdf -dEPSFitPage -dEPSCrop delme.eps
% epstopdf delme.eps

%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary directory
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('When finnished with the example, ');
disp('remove the temporary directory in which you are located!')
disp('>p=pwd;cd ..; rmdir(p,''s'');');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
