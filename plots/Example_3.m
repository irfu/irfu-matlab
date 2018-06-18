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
tint=[irf_time([2006 9 27 17 17 0]) irf_time([2006 9 27 17 24 0])]; % time interval

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1 % put to 0 if data already downloaded !!!!
    caa_download(tint,'C1_CP_FGM_5VPS');
	caa_download(tint,'CL_SP_AUX');
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
h=irf_plot(3,'newfigure'); % 3 subplots, remove middle (fast solution)
delete(h(2));				% remove middle panel
h(2)=[];					% remove handle
set(h(2),'position',get(h(2),'position')+[0.2 0.1 -0.2 0]); % move panel 2 up a bit

%%%%%%%%%%%%%%%%%%%%%%%
% top panel
hca=h(1);
% read data
B=irf_get_data('B_vec_xyz_gse__C1_CP_FGM_5VPS','caa','mat');
% plot
irf_plot(hca,B);
ylabel(hca,'B [nT] GSE');
irf_zoom(hca,'y',[-25 15])
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%
% top panel
hca=h(2);
% read data
B=irf_get_data('B_vec_xyz_gse__C1_CP_FGM_5VPS','caa','mat');
% plot
irf_plot(hca,B);
ylabel(hca,'B [nT] GSE');
irf_zoom(hca,'y',[-25 15])
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
irf_legend(hca,{'C1'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_pl_number_subplots(h);
irf_zoom(h(1),'x',tint);
% zoom figure 2 to smaller interval
tzoom=tint(1)+diff(tint)*[0.6 0.8];
irf_zoom(h(2),'x',tzoom);
% get satellite position
R=irf_get_data('sc_r_xyz_gse__CL_SP_AUX','caa','mat');
dR1=irf_get_data('sc_dr1_xyz_gse__CL_SP_AUX','caa','mat');
R1=irf_add(1,R,1,dR1);
irf_units % to get RE value
R1RE=irf_tappl(R1,'/Units.RE*Units.km'); % 
xx=get(gcf,'userdata');tst=xx.t_start_epoch;
xlab={'X (RE)','Y (RE)','Z (RE)'};
irf_timeaxis(h(1),tst,R1RE,xlab);
irf_timeaxis(h(1),'nodate');
% connect zoom in region
irf_plot_zoomin_lines_between_panels(h(1),h(2));
%
irf_legend(h(1),'Example 3',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% add interval mark
irf_pl_mark(h(1),tzoom,'blue')
irf_pl_mark(h(2),tzoom+[50 -30],'yellow')


%%%%%%%%%%%%%%%%%%%%%%%%
% to print the figure uncomment the lines below
%
% set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% print -dpng -painters Example_3.png;

%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary directory
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('When finnished with the example, ');
disp('remove the temporary directory in which you are located!')
disp('>p=pwd;cd ..; rmdir(p,''s'');');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!')


