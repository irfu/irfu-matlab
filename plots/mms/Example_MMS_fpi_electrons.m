%here it is shown how to plot several electron parameters
%NOTE: the fpi files available when this script was prepared are 10s
%resolution and NOT properly calibrated!
% by Sergio Toledo Redondo 2015/09/01

%% initialize database
mms.db_init('local_file_db','/data/mms');

%% select time interval
%tint  = irf.tint('2015-08-28T11:00:05Z/2015-08-28T18:00:00Z'); % define event time interval
tint  = irf.tint('2015-08-15T12:00:05Z/2015-08-15T15:00:00Z'); % define event time interval

ic=2;%spacecraft
%% Load data 

%if(exist('vex','var')==0);
    c_eval('B=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
    c_eval(' E=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
    %fpi data
    c_eval('ni=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint);',ic);
    c_eval('ne=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
    c_eval('vex=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_X_DSC'',tint);',ic);
    c_eval('vey=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Y_DSC'',tint);',ic);
    c_eval('vez=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Z_DSC'',tint);',ic);
    c_eval('tepar=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPara'',tint);',ic);
    c_eval('teper=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPerp'',tint);',ic);
    c_eval('Pexx=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_XX_DSC'',tint);',ic);
    c_eval('Peyy=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_YY_DSC'',tint);',ic);
    c_eval('Pezz=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_ZZ_DSC'',tint);',ic);
    c_eval('eEnSp_pX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pX'',tint);',ic);%positiveX
    c_eval('eEnSp_pY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pY'',tint);',ic);
    c_eval('eEnSp_pZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pZ'',tint);',ic);
    c_eval('eEnSp_mX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mX'',tint);',ic);%negativeX
    c_eval('eEnSp_mY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mY'',tint);',ic);
    c_eval('eEnSp_mZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mZ'',tint);',ic);
    c_eval('j_fpi_x=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_X_DSC'',tint);',ic);
    c_eval('j_fpi_y=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Y_DSC'',tint);',ic);
    c_eval('j_fpi_z=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Z_DSC'',tint);',ic);
%end

ve=irf.ts_vec_xyz(vex.time,[vex.data vey.data vez.data]);
Pev=irf.ts_vec_xyz(Pexx.time,[Pexx.data Peyy.data Pezz.data]);
j_fpi=irf.ts_vec_xyz(j_fpi_x.DEPEND_0.data,[j_fpi_x.data,j_fpi_y.data,j_fpi_z.data]);
t_e=irf.ts_vec_xy(tepar.DEPEND_0.data,[tepar.data teper.data]);

%% plot data



h = irf_plot(13,'newfigure');
h = irf_plot(13,'newfigure');
xSize=950; ySize=850;
set(gcf,'Position',[10 10 xSize ySize]); %position and size of the figure in the screen


hca = irf_panel('B'); %set(hca,'ColorOrder',mmsColors)
irf_plot(hca,B)
irf_legend(hca,{'B_x','B_y','B_z'},[0.02 0.1]);
ylabel(hca,[' B  ';'[nT]']);

hca = irf_panel('E');
irf_plot(hca,E)
ylabel(hca,['  E   ';'[mV/m]']);
irf_legend(hca,{'E_x','E_y','E_z'},[0.02 0.1]);
irf_zoom(hca,'y',[-20 20]);

hca = irf_panel('n'); 
irf_plot(hca,{ni,ne},'comp')
irf_legend(hca,{'n_i','n_e'},[0.02 0.1]);
ylabel(hca,['    n    '; '[cm^{-3}]'])

hca = irf_panel('ve');
irf_plot(hca,ve)
irf_legend(hca,{'v_{ex}','v_{ey}','v_{ez}'},[0.02 0.1]);
ylabel(hca,[' v_e  ';'[km/s]'])

 hca = irf_panel('te');
 irf_plot(hca,t_e);
 irf_legend(hca,{'T_{||}','T_{T}'},[0.02 0.1]);
 ylabel(hca,['Te ';'[K]']);

hca = irf_panel('Pe');
irf_plot(hca,Pev)
irf_legend(hca,{'P_{exx}','P_{eyy}','P_{ezz}'},[0.02 0.1]);
ylabel(hca,['     Pe     ';'[erg/cm^{3}]'])
    colormap('jet');

    
 hca = irf_panel('Jfpi'); 
irf_plot(hca,j_fpi)
irf_legend(hca,{'j_{x}','j_{y}','j_{z}'},[0.02 0.1]);
ylabel(hca,[' j  ';'[au]'])
    colormap('jet');   
    
    
%% 6 energy distributions    
hca=irf_panel('eEnSp_pX');
  specrec=struct('t',irf_time(eEnSp_pX.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_pX.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,+x}';'  [au]  ']);

hca=irf_panel('eEnSp_pY');
  specrec=struct('t',irf_time(eEnSp_pY.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_pY.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,+y}';'  [au]  ']);


hca=irf_panel('eEnSp_pZ');
  specrec=struct('t',irf_time(eEnSp_pZ.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_pZ.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,+z}';'  [au]  ']);

hca=irf_panel('eEnSp_mX');
  specrec=struct('t',irf_time(eEnSp_mX.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_mX.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,-x}';'  [au]  ']);

hca=irf_panel('eEnSp_mY');
  specrec=struct('t',irf_time(eEnSp_mY.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_mY.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,-y}';'  [au]  ']);


hca=irf_panel('eEnSp_mZ');

  specrec=struct('t',irf_time(eEnSp_mZ.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=eEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=eEnSp_mZ.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1E 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{e,-z}';'  [au]  ']);

%finish the plot
irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
%add_position(h(end),gsmR1)
xlabel(h(end),'')
title(h(1),['C' int2str(ic) ' ' tint.start.utc]);

