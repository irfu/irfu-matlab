%% Examples of MMS data usage
% Includes: 
% 1. Defining time series
% 2. Operations on time series
% 3. Plotting time series
% 4. MMS data structure
% 5. MMS position
% 6. Reading MMS E field data
% ...
%
% To open the example file and jump to some specific example execute
% > edit mmsdemo
% You should enable code folding in the editor!
%
% To execute demo run  
% > echodemo mmsdemo
%
%% Control log information
% How much information is shown can be controled by IRF.LOG settings.

% For running examples we disable showing additional log information

irf.log('off')

% Enable it afterwards by choosing log level, e.g. 
% > irf.log('critical') 

% Other levels are 'warning','notice','debug'
%% Example  1 Plot artificial data
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% showing exponentially growing wave and plot. It is good idead to get used 
% to using axis handles (variable 'h' in example). 

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T.tts - T.tts(1);                     % define relative time in s from start
x   = exp(0.001*(t)).*sin(2*pi*t/180);      % define function x(t)=exp(0.001(t-to))*sin(t-to)
TS1 = TSeries(T,x);                         % define TSeries object

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
%% Example  2 Plot multicomponent data
% Generate data with two components and plot in the same figure.
% Add legend text in lower left corner.

y   = exp(0.001*(t)).*cos(2*pi*t/180);	      % y(t)=exp(0.001(t-to))*cos(t-to)
TS2 = TSeries(T,[x y]);                       % define TSeries object

irf_plot(h,TS2)                               % plot in the same axis
irf_legend(h,{'X','Y'},[0.02 0.02])           % add legend text with the same colors as lines
%% Example  3 Work with data, zoom in plots. 
% Generate second data set that is a function of first.
% Plot both data sets in separate panels.
% Zoom in to smaller 30min time interval .

h    = irf_plot(2,'reset');		                 % empty figure with 2 panels
TS3 = TS2*1.2+2;		                           % TS2 = TS1*1.2 + 2	 

irf_plot(h(1),TS2);
irf_legend(h(1),{'X','Y'},[0.02 0.98],'fontsize',20)
irf_plot(h(2),TS3);
ylabel(h(2),'TS3 = TS2 * 1.2 +2 [nT]');

tintZoom = irf_time('2002-03-04T09:50:00Z/2002-03-04T10:01:00Z','utc>tint'); % 11 min interval
irf_zoom(h,'x',tintZoom);
irf_zoom(h,'y');
%% Example  4 Compare two data  
% Compare component-wise two datasets.
% Add one more label row with hours from the beginning of time interval

h=irf_plot({TS2,TS3},'comp');
ylabel(h(1),'B_X');
title(h(1),irf_time(TS2.t.tts(1),'tt>utc_yyyy-mm-dd'));
ylabel(h(2),'B_Y');
irf_legend(h(1),{'B2','Bnew=B2*1.2+2 '},[0.02 0.98],'fontsize',20)

TShours =TSeries(T,t/3600); % hours from beginning of time interval
irf_timeaxis(h(2),t(1),TShours,{'hours'})
irf_timeaxis(h(end),'nodate');
%% Example  5 Plot, different markers, mark intervals
% Plot using different markers or just some components of data. 
% irf_plot accepts all the parameters as normal matlab plot routine
% Mark some interesting time interval.

tint2=[irf_time([2008 03 01 10 11 0]) irf_time([2008 03 01 10 15 0])];
h=irf_plot(2,'reset');
irf_plot(h(1),TS2);
irf_plot(h(2),TS2,'r.','markersize',12);
ylabel(h(2),'B [nT] sc2');
irf_zoom(h,'y');
irf_pl_mark(h(2),tint2);
irf_legend(0,'Some additional info',[0,1],'color','r')
%% Example  6 Get Cluster Active Archive data
% This example requires that you are connected to internet!
% Get some real data from Cluster Active Archive. 
% First list what products are available having FGM_SPIN in name for given time
% interval
tmpDir = tempname;
mkdir(tmpDir);
cd(tmpDir);

tint=[irf_time([2003 02 04 18 0 0]) irf_time([2003 02 04 20 0 0])];
caa_download(tint,'list:*FGM_SPIN*')
%% Example  6 cont.
% Data exist from all 4 s/c, both in ISR2 and GSE reference frame.
% Download C1_CP_FGM_SPIN that has B in GSE reference frame.
% Data are saved as cdf files in subdirectories under ./CAA directory.
% You can see what data have been downloaded from CAA and load those into memory
% using caa_load command. When data are loaded into the memory they are in dataobj
% form that includes all the metadata information within cdf file. Just enter
% the databoj name to see the variables within it. To get other meta information
% you can try 
% > C1_CP_FGM_SPIN.GlobalAttributes
% > C1_CP_FGM_SPIN.VariableAttributes

caa_download(tint,'C1_CP_FGM_SPIN');	
caa_load list;
caa_load C1_CP_FGM_SPIN
C1_CP_FGM_SPIN
%% Example  6 cont.
% If you just want to see some variable you can directly plot it knowing its
% name. If you want to manipulate variables, you will need to extract them from
% dataobject. In example they are extracted in matlab format, first column time
% and other components are vector components.
% If vector magnitude is missing, you can fast calculate it and add as last
% comoponent using irf_abs. 

irf_plot('B_vec_xyz_gse__C1_CP_FGM_SPIN');
B1 = c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_SPIN','mat');
B1 = irf_abs(B1);
h=irf_plot(1,'newfigure');
irf_plot(h,B1);
irf_legend(h,{'Bx','By','Bz','B'},[0.02 0.98])
%% Example  7 Cluster location
% Cluster location. The position is downloaded from internet, isdat server in
% Uppsala. If you want to use position also offline, you have to download it
% from CAA. Command below would download position files for all 4 s/c
%
% > caa_download(tint,'C?_CP_AUX_POSGSE_1M');
%
% Please note that under Menu 'Options' there are different ways how
% the satellite configuration can be visualized. 

c_pl_sc_conf_xyz
%% Example  8 Both B-field and particle data (spectra)
% Example that plot both magnetic field and ions spectrogram. 
% First download ion data. You can list different data using
% > caa_download(tint,'list:*CIS*')

caa_download(tint,'*CODIF_O1_1D*PEF')

h=irf_plot(2,'newfigure');
irf_plot(h(1),'B_vec_xyz_gse__C1_CP_FGM_SPIN');
irf_plot(h(2),'flux__C1_CP_CIS_CODIF_O1_1D_PEF')

set(h(2),'yscale','log')
irf_zoom(h,'x',tint)
irf_plot_axis_align(h)		% align X-limits of axes boxes 
%% Example  8 cont. (particle moments)
% Download also moments. Use different colorscale for spectra. 
caa_download(tint,'*CODIF*O*MOMENTS*')

clf;
h=irf_plot(3);
irf_plot(h(1),'B_vec_xyz_gse__C1_CP_FGM_SPIN');
irf_plot(h(2),'flux__C1_CP_CIS_CODIF_O1_1D_PEF')
irf_colormap(h(2),'space')
irf_plot(h(3),'velocity__C1_CP_CIS_CODIF_HS_O1_MOMENTS')
irf_zoom(h(3),'y',[-300 300]);
set(h(2),'yscale','log')
irf_zoom(h,'x',tint);
irf_plot_axis_align(h)

irf_legend(h(3),{'VX','VY','VZ'},[0.02 0.05])
%% Example 10 MVAR Minimum Variance Analysis 
% Get the data from file and run minimum variance analysis. 
% You select interval by clicking with mouse (red selection shows whether you
% have to select the beginning or end of time interval). 
% Calculate MVA by pressing menu "Recalculate"
%
% Data files are downloaded from CAA in Example 6

B1 = c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_SPIN','mat');
B1 = irf_abs(B1);
irf_minvar_gui(B1)
%% Example 11 Multispacecraft timing analysis
% Run interactive timing analysis that allows you to estimate phase velocity of
% boundaries from time of boundary crossings. When working with all 4 spacecraft
% data very useful is routine c_eval, that allows to execute commands on list of
% satellites (default 1..4), see the help. 
% You will need to be connected to internet to get sc position, unless you have
% downloaded the position data before as described in Example 7. 
%
% Once window is open you can click menu "Click times" and by pointer mark the
% boundaries that you want to time. You can express the result in distance
% clicking menu "Distance". Reference satellite is not time shifted.
%
% You can also enter manually either time offsets or boundary speeds.
%

% Download data for example magnetopause crossing
tint=[irf_time([2002 03 14 00 20 0]) irf_time([2002 03 14 00 25 0])];
caa_download(tint,'C*_CP_FGM_SPIN','overwrite');	
% C* to download all sc data, 'overwrite' to remove previous data

c_eval('B?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_SPIN'',''mat'',''file'');');
c_eval('B?=irf_abs(B?);');
c_4_v_gui('B?')

