% A routine to compute and plot reduced electron distributions from FPI
%
% Compare with: Wilder, F. D., et al. (2016), GRL, 43, 5909?5917,
% doi:10.1002/2016GL069473.
%
% The Example is fairly slow. Approx 2 min.
%
% Originally written for 1SC at two time intervals by A. Johlander
% Adapted for 2SC (or all SC) at a single time by J. White


%% Set parameters and get data
tStart = tic;
% time interval
tint = irf.tint('2017-07-06T00:54:12.000/2017-07-06T00:54:20.00Z');
% time to plot reduced electron distribution
t1 = irf.time_array('2017-07-06T00:54:14.000');

% sc numbers to compare in time series
sc = [1,4];

% if you want to compare all SC reduced distributions, allSC = true
% else compare SC in sc (above)
allSC = true;

% plot either 'reduced', 'omni', or 'pitchangle' distributions
% plotDistType = 'reduced'; % could have cases or make this arg to funct

% pitch angle range to plot
pitchRange = [0,30];

% color/y-limit
clim = 10.^[-5,-2.5]; % s m^-4

% define velocity grid
vg = linspace(-60e3,60e3,100); % km/s

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 1e3;

% velocity limit in plot
vlim = [min(vg),max(vg)]; % km/s

% matlab colours
colours = [0.0000    0.4470    0.7410;...
           0.8500    0.3250    0.0980];
allSC_colours = [0   0   0;...
                 0.8 0   0;...
                 0   0.8 0;...
                 0   0   0.8];

% get distribution functions
c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint));',1:4)
c_eval('ePDist? = ePDist?.tlim(tint);',1:4)

% get magnetic field in DMPA (since ePDist is in DMPA)
c_eval('B? = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);',1:4)

% Get the spacecraft potential (in this example, it is not so important)
c_eval('scPot? = mms.get_data(''V_edp_brst_l2'',tint,?);',1:4)

% remove flux from bottom two energy levels to make it more like in the
% paper
c_eval('ePDist?.data(:,1:2,:,:) = 0;',1:4)

% get index for the time, number of times to plot
c_eval('it = interp1(ePDist?.time.epochUnix,1:length(ePDist?),t1.epochUnix,''nearest'');',sc(1))
numTimes = 5; % use an odd number
halfTime = round((numTimes-1)/2);

%% Reduce distribution
tic
% reduced distribution along B
c_eval('f1D_? = ePDist?.reduce(''1D'',B?,''vg'',vg,''nMC'',nMC,''scpot'',scPot?);',1:4)
toc

%% Plot reduced distribution as a time series
% make figure
h = irf_plot(5,'newfigure');

% plot magnetic fields for sc
for k=[1,2]
  hca = irf_panel(h,['Bxyz ',num2str(k)]);
  c_eval('irf_plot(hca,B?);',sc(k))
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,['B_{dmpa,',num2str(sc(k)),'} [nT]'])
  irf_legend(hca,{'Bx';'By';'Bz'},[1.02,0.9])
end

% Plot reduced distributions
for k=[1,2]
  hca = irf_panel(h,['pdist ',num2str(k)]);
  c_eval('[~,hcb] = irf_spectrogram(hca,f1D_?.specrec(''1D_velocity''));',sc(k))
  hcb.Label.String = 'log_{10}f^{1D}_e [s m^{-4}]';
  ylabel(hca,['V_{||,',num2str(sc(k)),'} [km/s]'])
  colormap('jet')
  irf_zoom(hca,'y',[min(vg),max(vg)])
end

% subtracting distributions from each other
c_eval('f1D_?! = f1D_?; f1D_?!.data = f1D_?.data - f1D_!.data;',sc(1),sc(2))

hca = irf_panel(h,'pdist diff');
[~,hcb] = irf_spectrogram(hca,f1D_14.specrec('1D_velocity'), 'lin');
hcb.Label.String = 'log_{10}F_e [s m^{-4}]';
ylabel(hca,['\DeltaV_{||,',num2str(sc(1)),',',num2str(sc(2)),'} [km/s]'])
colormap('jet')
irf_zoom(hca,'y',[min(vg),max(vg)])
hca.CLim = [-10^(-3),10^(-3)];

% plot vertical black line to show when distribution slice is
vertLine = irf_pl_mark(h,t1,'k');

% align axes
irf_plot_axis_align(h)
irf_zoom(h,'x',tint)
h(1).Title.String = ['MMS ',num2str(sc(1)),' & ',num2str(sc(2)),' Data'];

%% Plot reduced distributions as a function of parallel velocity

% initiate figure
hca = irf_plot(1,'newfigure');
hold(hca,'on')

if ~allSC % plot only SC in sc array
  if halfTime > 0
    % plot multiple individual times
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?.data(it+[-!:1:!],:),'':'',''Color'',colours(2,:),''linewidth'',1)',sc(1),halfTime)
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?.data(it+[-!:1:!],:),'':'',''Color'',colours(1,:),''linewidth'',1)',sc(2),halfTime)

    % plot time-averaged reduced distribution function (i.e. integrated over
    % time and divided by scalar factor)
    c_eval('plot(hca,f1D_?(it).depend{1},mean(f1D_?.data(it+[-!:1:!],:)),''Color'',colours(2,:),''linewidth'',2)',sc(1),halfTime)
    c_eval('plot(hca,f1D_?(it).depend{1},mean(f1D_?.data(it+[-!:1:!],:)),''Color'',colours(1,:),''linewidth'',2)',sc(2),halfTime)

    % create string to append to title string
    titleStringAdd = [' (time-averaged for ', num2str(numTimes), ' times)'];
  else
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?(it).data,''Color'',colours(2,:),''linewidth'',2)',sc(1))
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?(it).data,''Color'',colours(1,:),''linewidth'',2)',sc(2))

    % create string to append to title string
    titleStringAdd = ' (not time-averaged)';
  end

  % legend shows SC
  hca.ColorOrder = flipud(colours(1:2,:)); % set colour order for legends
  irf_legend(hca,{['MMS',num2str(sc(1))];['MMS',num2str(sc(2))];},[0.98,0.98])

else % plot all SC
  if halfTime > 0
    % plot multiple individual times as dashed lines
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?.data(it+[-!:1:!],:),'':'',''Color'',allSC_colours(?,:),''linewidth'',1)',1:4,halfTime)

    % plot time-averaged reduced distribution function (i.e. integrated over
    % time and divided by scalar factor) as bold lines
    c_eval('plot(hca,f1D_?(it).depend{1},mean(f1D_?.data(it+[-!:1:!],:)),''Color'',allSC_colours(?,:),''linewidth'',2)',1:4,halfTime)

    % create string to append to title string
    titleStringAdd = [' (time-averaged for ', num2str(numTimes), ' times)'];
  else
    c_eval('plot(hca,f1D_?(it).depend{1},f1D_?(it).data,''Color'',allSC_colours(?,:),''linewidth'',2)',1:4)

    % create string to append to title string
    titleStringAdd = ' (not time-averaged)';
  end

  % legend shows SC
  hca.ColorOrder = allSC_colours; % set colour order for legend
  irf_legend(hca,{'MMS1';'MMS2';'MMS3';'MMS4'},[0.98,0.98])

end

% scales
hca.YScale = 'log';
hca.YLim = clim;
hca.XLim = vlim;
% labels
xlabel(hca,'V_{||} [km/s]')
ylabel(hca,'f^{1D}_e [s m^-^4]')
c_eval('timeString = ePDist?(it).time.toUtc;',sc(1))
timeOfDay = timeString(12:23); % get time string substring 'HH:MM:SS.mmm'
hca.Title.String = [timeOfDay, titleStringAdd];

%% Plot omni distributions as a function of energy

% initiate figure
hca = irf_plot(1,'newfigure');
hold(hca,'on')

if ~allSC % plot only SC in sc array
  if halfTime > 0
    % plot multiple individual times
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?.omni.data(it+[-!:1:!],:),'':'',''Color'',colours(2,:),''linewidth'',1)',sc(1),halfTime)
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?.omni.data(it+[-!:1:!],:),'':'',''Color'',colours(1,:),''linewidth'',1)',sc(2),halfTime)

    % plot time-averaged reduced distribution function (i.e. integrated over
    % time and divided by scalar factor)
    c_eval('plot(hca,ePDist?(it).depend{1},mean(ePDist?.omni.data(it+[-!:1:!],:)),''Color'',colours(2,:),''linewidth'',2)',sc(1),halfTime)
    c_eval('plot(hca,ePDist?(it).depend{1},mean(ePDist?.omni.data(it+[-!:1:!],:)),''Color'',colours(1,:),''linewidth'',2)',sc(2),halfTime)

    % create string to append to title string
    titleStringAdd = [' (time-averaged for ', num2str(numTimes), ' times)'];
  else
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?(it).omni.data,''Color'',colours(2,:),''linewidth'',2)',sc(1))
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?(it).omni.data,''Color'',colours(1,:),''linewidth'',2)',sc(2))

    % create string to append to title string
    titleStringAdd = ' (not time-averaged)';
  end

  % legend shows SC
  hca.ColorOrder = flipud(colours(1:2,:)); % set colour order for legends
  irf_legend(hca,{['MMS',num2str(sc(1))];['MMS',num2str(sc(2))];},[0.98,0.98])

else % plot all SC
  if halfTime > 0
    % plot multiple individual times as dashed lines
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?.omni.data(it+[-!:1:!],:),'':'',''Color'',allSC_colours(?,:),''linewidth'',1)',1:4,halfTime)

    % plot time-averaged reduced distribution function (i.e. integrated over
    % time and divided by scalar factor) as bold lines
    c_eval('plot(hca,ePDist?(it).depend{1},mean(ePDist?.omni.data(it+[-!:1:!],:)),''Color'',allSC_colours(?,:),''linewidth'',2)',1:4,halfTime)

    % create string to append to title string
    titleStringAdd = [' (time-averaged for ', num2str(numTimes), ' times)'];
  else
    c_eval('plot(hca,ePDist?(it).depend{1},ePDist?(it).omni.data,''Color'',allSC_colours(?,:),''linewidth'',2)',1:4)

    % create string to append to title string
    titleStringAdd = ' (not time-averaged)';
  end

  % legend shows SC
  hca.ColorOrder = allSC_colours; % set colour order for legend
  irf_legend(hca,{'MMS1';'MMS2';'MMS3';'MMS4'},[0.98,0.98])

end

% scales
hca.YScale = 'log';
hca.XScale = 'log';
%hca.YLim = clim;
%hca.XLim = vlim;
% labels
xlabel(hca,'E_e [eV]')
ylabel(hca,'f_e [s^3 m^{-6}]')
c_eval('timeString = ePDist?(it).time.toUtc;',sc(1))
timeOfDay = timeString(12:23); % get time string substring 'HH:MM:SS.mmm'
hca.Title.String = [timeOfDay, titleStringAdd];

%% When running full script, display time taken
tStop = toc(tStart);
disp(['Time taken to run example ~', ...
       num2str(round(tStop)), ...
       ' seconds.'])
