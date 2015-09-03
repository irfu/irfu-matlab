%here it is shown how to plot several ion parameters
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
%file='/DATA/mms/mms?/fpi/fast/sitl/2015/08/mms?_fpi_fast_sitl_20150815115000_v0.0.0.cdf';

%if(exist('B','var')==0) %to run it only once
    c_eval('B=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
    c_eval(' E=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
    %fpi data
    c_eval('ni=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint);',ic);
    c_eval('ne=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
    c_eval('vix=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_X_DSC'',tint);',ic);
    c_eval('viy=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Y_DSC'',tint);',ic);
    c_eval('viz=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Z_DSC'',tint);',ic);
    c_eval('tipar=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_DIStempPara'',tint);',ic);
    c_eval('tiper=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_DIStempPerp'',tint);',ic);
    c_eval('Pixx=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISpress_XX_DSC'',tint);',ic);
    c_eval('Piyy=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISpress_YY_DSC'',tint);',ic);
    c_eval('Pizz=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISpress_ZZ_DSC'',tint);',ic);
    c_eval('iEnSp_pX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pX'',tint);',ic);%positiveX
    c_eval('iEnSp_pY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pY'',tint);',ic);
    c_eval('iEnSp_pZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pZ'',tint);',ic);
    c_eval('iEnSp_mX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mX'',tint);',ic);%negativeX
    c_eval('iEnSp_mY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mY'',tint);',ic);
    c_eval('iEnSp_mZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mZ'',tint);',ic);
    c_eval('j_fpi_x=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_X_DSC'',tint);',ic);
    c_eval('j_fpi_y=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Y_DSC'',tint);',ic);
    c_eval('j_fpi_z=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Z_DSC'',tint);',ic);
%end


%fpi data handling
vi=irf.ts_vec_xyz(vix.time,[vix.data viy.data viz.data]);
Piv=irf.ts_vec_xyz(Pixx.time,[Pixx.data Piyy.data Pizz.data]);
j_fpi=irf.ts_vec_xyz(j_fpi_x.DEPEND_0.data,[j_fpi_x.data,j_fpi_y.data,j_fpi_z.data]);
t_i=irf.ts_vec_xy(tipar.DEPEND_0.data,[tipar.data tiper.data]);




%% plot data


h = irf_plot(13,'newfigure');
xSize=950; ySize=850;
set(gcf,'Position',[10 10 xSize ySize]); %position and size of the figure in the screen

hca = irf_panel('B');
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

hca = irf_panel('vi'); 
irf_plot(hca,vi)
irf_legend(hca,{'v_{ix}','v_{iy}','v_{iz}'},[0.02 0.1]);
ylabel(hca,[' v_i  ';'[km/s]'])

 hca = irf_panel('ti'); 
 irf_plot(hca,t_i);
 irf_legend(hca,{'T_{||}','T_{T}'},[0.02 0.1]);
 ylabel(hca,['Ti ';'[K]']);

hca = irf_panel('Pi'); 
irf_plot(hca,Piv)
irf_legend(hca,{'P_{ixx}','P_{iyy}','P_{izz}'},[0.02 0.1]);
ylabel(hca,['     Pi     ';'[erg/cm^{3}]'])
colormap('jet');

    
hca = irf_panel('Jfpi');
irf_plot(hca,j_fpi)
irf_legend(hca,{'j_{x}','j_{y}','j_{z}'},[0.02 0.1]);
ylabel(hca,[' j  ';'[au]'])
colormap('jet');   
    
    
%% 6 energy distributions    
hca=irf_panel('iEnSp_pX');
  specrec=struct('t',irf_time(iEnSp_pX.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_pX.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,+x}';'  [au]  ']);

hca=irf_panel('iEnSp_pY');
  specrec=struct('t',irf_time(iEnSp_pY.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_pY.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,+y}';'  [au]  ']);


hca=irf_panel('iEnSp_pZ');
  specrec=struct('t',irf_time(iEnSp_pZ.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_pZ.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,+z}';'  [au]  ']);

hca=irf_panel('iEnSp_mX');
  specrec=struct('t',irf_time(iEnSp_mX.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_mX.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,-x}';'  [au]  ']);

hca=irf_panel('iEnSp_mY');
  specrec=struct('t',irf_time(iEnSp_mY.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_mY.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,-y}';'  [au]  ']);


hca=irf_panel('iEnSp_mZ');

  specrec=struct('t',irf_time(iEnSp_mZ.DEPEND_0.data,'ttns>epoch'));
    %specrec.f=iEnSp_pX.DEPEND_1.data;%energy levels
    specrec.f=0:31;%energy levels
    specrec.p=iEnSp_mZ.data(:,1:32);%data matrix
    specrec.f_label='a';
    specrec.p_label='Energy';
    irf_spectrogram(hca,specrec);
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0, 5])
ylabel(hca,['E_{i,-z}';'  [au]  ']);

%finish the plot
irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
%add_position(h(end),gsmR1)
xlabel(h(end),'')
title(h(1),['C' int2str(ic) ' ' tint.start.utc]);







