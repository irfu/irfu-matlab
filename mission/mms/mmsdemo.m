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
TS1 = irf.ts_scalar(T,x);                   % define scalar TSeries object

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
%% Example  2 Plot multicomponent data
% Generate data with two components and plot in the same figure.
% Add legend text in lower left corner.

y   = exp(0.001*(t)).*cos(2*pi*t/180);	      % y(t)=exp(0.001(t-to))*cos(t-to)
TS2 = irf.ts_vec_xy(T,[x y]);                 % define vector TSeries object

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

%tintZoom = irf_time('2002-03-04T09:50:00Z/2002-03-04T10:01:00Z','utc>tint'); % 11 min interval
tintZoom = irf.tint('2002-03-04T09:50:00Z/2002-03-04T10:01:00Z'); % 11 min interval
irf_zoom(h,'x',tintZoom);
irf_zoom(h,'y');
%% Example  4 Compare two data  
% Compare component-wise two datasets.
% Add one more label row with hours from the beginning of time interval

h=irf_plot({TS2,TS3},'comp');
ylabel(h(1),'B_X');
title(h(1),TS2.time(1).utc('yyyy-mm-dd'));
ylabel(h(2),'B_Y');
irf_legend(h(1),{'B2','Bnew=B2*1.2+2 '},[0.02 0.98],'fontsize',20)

TShours =irf.ts_scalar(T,t/3600); % hours from beginning of time interval
irf_timeaxis(h(2),t(1),TShours,{'hours'})
irf_timeaxis(h(end),'nodate');
%% Example  5 Plot, different markers, mark intervals
% Plot using different markers or just some components of data. 
% irf_plot accepts all the parameters as normal matlab plot routine
% Mark some interesting time interval.

tint2=[irf_time([2002 03 04 9 50 0]) irf_time([2002 03 04 9 55 0])];
h=irf_plot(2,'reset');
irf_plot(h(1),TS2);
irf_plot(h(2),TS2,'r.','markersize',12);
ylabel(h(2),'B [nT] sc2');
irf_zoom(h,'y');
irf_pl_mark(h(2),tint2);
irf_legend(0,'Some additional info',[0,1],'color','r')
%% Generate vector times series
TS2 = irf.ts_vec_xyz(T,[x y y]);             % define TSeries object
irf_minvar_gui(TS2)

%% MMS data structure

%% Example  6 Get MMS data

%% Example  7 MMS location
% MMS location. 

%mms_configuration

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

%c_4_v_gui('B?')
