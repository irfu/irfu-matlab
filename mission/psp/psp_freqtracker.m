%                               INFORMATION                               %
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
---------------------------------------------------------------------------
                                DESCRIPTION
---------------------------------------------------------------------------
A simple funcion following a specific frequency based on the increased amplitude
(in this example the plasma line) at the nearby bins spectrogram (e.g. from PSP lfr) 

Below you can find a detailed description of input/output and possible
future updates/implementations.

REQUIREMENTS: irfu-matlab (pspdevel), setargs, set up datastore
path for psp_load

%All outputs are optional.
%All inpuits are optional but you need to provide either a day or a start
%and stop period.
e.g. 
1) start = [2020 06 09 03 00 00] ; stop =  [2020 06 09 04 00 00];
2) fullday = [2020 06 09] ;
---------------------------------------------------------------------------
                                INPUT
---------------------------------------------------------------------------
'initialfreq'           - bin for initial frequency (default = 30)
'minfreq'               - minimum frequency bin (default = initialfreq-10)
'maxfreq'               - minimum frequency bin (default = initialfreq+10)
'binrange'              - bin range to search for higher amplitude (default = 2)
'tolerancelevel'        - rough estimation of background level to avoid changing
                          too frequenctly. See code below on how it is computed.
'recoverer'             - assisting parameter to return to a logical range of
                          values in case of "spikes" (default = 5)
'generateplot'          - generate plot for spectrum (default = false)
'metricsdesplay'        - display metrics of routine (default = true)
'initialdensitydata'    - assist initial frequency finder with spi data
                          (default = true
'densitycomparisonplot' - generate density comparison plots (default =
                           alse)

---------------------------------------------------------------------------
                                OUTPUT
---------------------------------------------------------------------------
[A,B,C,D]

A = Structure with plasma line TS and weighted plasma line TS
B = Strcture with density TS and weighted density TS
C = Structure with low resolution frequency & amplitude from lfr
D = Timeseries of density from spi data

---------------------------------------------------------------------------
                                TODO
---------------------------------------------------------------------------

TOP PRIORITY

- Write code that works with hires data from rfs (ONGOING)

HIGH PRIORITY

- Establishing background tolerance level changing over time in case
of different fluxes or density in close proximity.

- Re-visit algorithm to optimize computational time
and add a few callbacks in case of errors.

LOW PRIORITY

- Quality index based on differences for other instruments. (in other
words, when plasma line is hard to read, trust the other instruments)

- Allow dynamic hyperparameter optimization within the for loop (changing
the default parameters when things are changing drastically)

- Extra outputs/inputs (saving plots, etc.)

- Update with current updated irfu_matlab routines

---------------------------------------------------------------------------
                                EXAMPLES 
                          see irfu/plots/psp/
---------------------------------------------------------------------------

start = [2020 06 09 03 00 00] ;
stop =  [2020 06 09 04 00 00];
fullday = [2020 06 09] ;

plasmaline = FreqTracking(start,stop);
[plasmaline,DensityTS,pspdata,InitialDensityTS] = psp_freqtracker(start,stop,'initialdensitydata',true,'generateplot',true,'densitycomparisonplot',true);
[plasmaline,~,pspdata,~] = psp_freqtracker(fullday,'generateplot',false);

---------------------------------------------------------------------------

-Savvas Raptis & Sabrina Tigik Ferrao
KTH, Royal Institute of Technology
savvra@kth.se
2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

function [varargout] = psp_freqtracker(varargin)

%% Add units from irf
Units = irf_units; % read in standard units
Me = Units.me;
Mp = Units.mp;
c = Units.c;
e = Units.e;
epso = Units.eps0;
mu0 = Units.mu0;
Mp_Me = Mp/Me; % ratio of proton and electron mass;


%% Organizing input parameters and desirable output
if nargin == 0
    error("Provide date or start and stop time")
end

if  nargin > 0
    if length(varargin{1}) < 4 
        defaultstarttime = varargin{1};
        defaultstarttime = [defaultstarttime,00];
        defaultstoptime = varargin{1};
        defaultstoptime = [defaultstoptime,00];
        varargin = ['fullday' varargin(1) 'start' varargin(1) 'stop' varargin(1:end)];
    elseif length(varargin{1}) > 3 && nargin>1
        defaultstarttime = varargin{1};
        defaultstoptime =  varargin{2};
        varargin = ['fullday' varargin(1) 'start' varargin(1) 'stop' varargin(2:end)];
    else
        error("Provide date or start and stop time")
    end
end


SpectrumData = psp_load([],'pl_rfs_lfr', [defaultstarttime(1) defaultstarttime(2) defaultstarttime(3)], [defaultstoptime(1) defaultstoptime(2) defaultstoptime(3)]); %
amplTS = SpectrumData{1};
freqTS = SpectrumData{2};
amplTSHR = SpectrumData{3}; %todo
freqTSHR = SpectrumData{4}; %todo

% Establishing some default parameters that worked on the example case in case
% different ones are not available
defaultrecoverer = 5; %Searching for maximum or minimum frequency in case of "spike"
defaultinitialfreq = 30; % Some starting point, will search for better position if initialdensitydata is given
defaulttolerancelevel_weight = 1/2; % weight can range from (0,1) otherwise results will be very noisy
defaulttolerancelevel = defaulttolerancelevel_weight * min(amplTS.data(1,defaultinitialfreq)-amplTS.data(1,[defaultinitialfreq-2:defaultinitialfreq-1 defaultinitialfreq+1:defaultinitialfreq+2]));
defaultbinrange = 2; % bin range to search for higher amplitude

defaultargs = {'fullday',defaultstarttime,'start', defaultstarttime, 'stop', defaultstoptime,...
    'initialfreq', defaultinitialfreq, 'minfreq',defaultinitialfreq-10,'maxfreq',defaultinitialfreq+10, ...
    'binrange',defaultbinrange, 'tolerancelevel',defaulttolerancelevel','recoverer',defaultrecoverer,'generateplot',false, ...
    'metricsdesplay',true, 'initialdensitydata',true, 'densitycomparisonplot',false};


disp('The arguments choosen are:') 
argstruct   = setargs(defaultargs, varargin) % showing arguments used


if size(start,2) <4  %% Small hotfix to allow different initial inputs
    start = [start 00] ;
    stop = [stop 24] ;
end 

start = [start(1) start(2) start(3) start(4) 00 00] ; %Small hotfix on times to avoid a bug with minutes being different than zero
stop =  [stop(1) stop(2) stop(3) stop(4) 00 00]; %Small hotfix on times to avoid a bug with minutes being different than zero
startepochtt = irf_time(start,'vector>epochtt');
stopepochtt = irf_time(stop,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

testrfs_lfr_v1v2_freq = tlim(freqTS,tint); %Splitting the dataset
testrfs_lfr_v1v2 = tlim(amplTS,tint);  %Splitting the dataset

rfs_lfr_v1v2_spec=struct('t', testrfs_lfr_v1v2_freq.time.epochUnix);
rfs_lfr_v1v2_spec.p = testrfs_lfr_v1v2.data(:,1:size(testrfs_lfr_v1v2.data,2));
rfs_lfr_v1v2_spec.f = single(testrfs_lfr_v1v2_freq.data(:,1:size(testrfs_lfr_v1v2_freq.data,2)));
rfs_lfr_v1v2_spec.p_label={['log10 (' testrfs_lfr_v1v2.units ')']};
rfs_lfr_v1v2_spec.f_label={testrfs_lfr_v1v2_freq.units};

if initialdensitydata == true
    AssistingData = psp_load([],'spi', [start(1) start(2) start(3)], [start(1) start(2) start(3)]); % Download 1 day of data to local datastore path (edit datastore to change or add first argument with a path)    
    SubAssistingData =  tlim(AssistingData{2},tint);
    defaultinitialfreqvalue =  sqrt((SubAssistingData.data(1)*1e6)*e^2/Me/epso)/2/pi; % 9e3*sqrt(initialdensitydata.data(1));
    PlasmaLineFromFrequency =  sqrt((SubAssistingData*1e6)*e^2/Me/epso)/2/pi;
    
    [~,initialfreq] = min(abs(freqTS.data(1,:) - defaultinitialfreqvalue));
end

%--------------------------------------------------------------------------
%% Actual function starts below.
%--------------------------------------------------------------------------


%% Assisting arrays for the tracker below
DataFrequency=[];
DataFrequencyWeightedAverage=[];
frequencytracker=[];
TrackerNumber=0;
spikecounter=0;

%% Tracker of frequency

for i=1:length(testrfs_lfr_v1v2.time)
    
    MaximumAmplitude = testrfs_lfr_v1v2.data(i,initialfreq-binrange:initialfreq+binrange);
    
    [~, maxindex] = find(max(MaximumAmplitude)==MaximumAmplitude);
    
    if maxindex ~= binrange+1 && (max(MaximumAmplitude) - testrfs_lfr_v1v2.data(i,initialfreq))> tolerancelevel
        %disp(['changing initial frequency from ', num2str(initialfreq), ' to ', num2str(initialfreq-binrange-1+maxindex)])
        initialfreq = initialfreq-binrange-1+maxindex;
        frequencytracker = [frequencytracker, initialfreq ] ;
        TrackerNumber = TrackerNumber +1;
    end
    
    DataFrequency = [DataFrequency,testrfs_lfr_v1v2_freq.data(i,initialfreq)];
    NormalizedWeightedAmplitudeArray = normalize(testrfs_lfr_v1v2.data(i,initialfreq-1:initialfreq+1));
    WeightedResult = sum(testrfs_lfr_v1v2_freq.data(i,initialfreq-1:initialfreq+1).*NormalizedWeightedAmplitudeArray);
    DataFrequencyWeightedAverage = [DataFrequencyWeightedAverage,WeightedResult];
    
    if initialfreq < minfreq || initialfreq > maxfreq
        spikecounter = spikecounter + 1;
        warning('Possible spike encounter: Re-evaluate the tracking')
        tempinitialfreq = initialfreq;
        if initialfreq > maxfreq
            initialfreq = min(frequencytracker(TrackerNumber-recoverer:end));
        else initialfreq = max(frequencytracker(TrackerNumber-recoverer:end));
        end
        disp(['changing initial frequency from ', num2str(tempinitialfreq) ' to ',num2str(initialfreq)])
    end
end

%% Displaying metrics

if metricsdesplay == true
    disp('Metrics of proccedure:')
    disp('----------------------')
    disp(['Tracking frequency for: ' 'TODO'])
    disp(['Small variations of initial frequency: ' num2str(TrackerNumber)])
    disp(['Large variations / Spikes of initial frequency: ' num2str(spikecounter)])
    if spikecounter>5
        disp('TIP: If there are a lot of spikes either try a shorter period or change the tolerance/bins');
    end
end
%% Plotting

PlasmaLine = TSeries(testrfs_lfr_v1v2_freq.time,DataFrequency');
WeightedPlasmaLine = TSeries(testrfs_lfr_v1v2_freq.time,DataFrequencyWeightedAverage');

DensityData = ((((PlasmaLine.data*2*pi).^2)*Me*epso)/e^2)/1e6;
WeightedDensityData = ((((WeightedPlasmaLine.data*2*pi).^2)*Me*epso)/e^2)/1e6;
DensityTS = TSeries(testrfs_lfr_v1v2_freq.time,DensityData); %density timeseries
WeightedDensityTS = TSeries(testrfs_lfr_v1v2_freq.time,WeightedDensityData); %density timeseries
WeightedDensityDataSmoothed = movmean(WeightedDensityData,3);
WeightedDensityTSSmoothed =  TSeries(testrfs_lfr_v1v2_freq.time,WeightedDensityDataSmoothed); %density timeseries

if generateplot == true
    h = irf_plot(1,'newfigure');
    [hcc, hcd] = irf_spectrogram(h(1),rfs_lfr_v1v2_spec,'donotfitcolorbarlabel');
    set(h(1),'Yscale', 'log')
    irf_zoom(h(1),'y',[7e4 1.5e6])
    hcc.YScale = 'log';
    hcc.YTick = [1e5 2e5 2.5e5 1e6];
    hcc.ColorScale = 'log';
    set(hcd,'TickDir','out')
    caxis(h(1),[-16, -13]);
    colormap(h(1),'jet');
    hold(hcc,'on')
    %hb = irf_plot(h(1),PlasmaLine,'color','white','LineWidth',1.5);
    hc = irf_plot(h(1),WeightedPlasmaLine,'--','color','white','LineWidth',1.5);
    if initialdensitydata == true
        hd = irf_plot(h(1),PlasmaLineFromFrequency,'color','black','LineWidth',1.5);
    end
    l1 = ylabel(hcc,'Hz','FontSize',11);
    legend(hcc,{'lfr v1v2','Weighted plasma line','Plasma Line SPI'},'textcolor','white','color','none','Interpreter','tex','FontSize',15,'location','northeast','Box','off')
    hold(hcc,'off')
end

if densitycomparisonplot == true
    k = irf_plot(1,'newfigure');
    hca = irf_panel('density');
    hold(hca,'on');
    irf_plot(WeightedDensityTSSmoothed,'LineWidth',1.5)
    if initialdensitydata == true
        irf_plot(SubAssistingData)
        ylim([min([DensityTS.data])-50,max([SubAssistingData.data])+50])
    else
        ylim([min([DensityTS.data])-50,max([DensityTS.data])+50])
    end
    irf_legend(hca,{'n_{weighted,smoothed}','n_{SPI}'},[0.95 0.97],'FontSize',11);
    l1 = ylabel(hca,'n_{e} (cm^-3)','FontSize',11);
    set(l1,'interpreter','tex');
    hold(hca,'off');
end


if nargout  == 0
    return
elseif nargout == 1
    varargout {1} = {PlasmaLine,WeightedPlasmaLine};
elseif nargout == 2
    varargout {1} = {PlasmaLine,WeightedPlasmaLine};
    varargout {2} = {DensityTS,WeightedDensityTSSmoothed};
elseif nargout == 3
    varargout {1} = {PlasmaLine,WeightedPlasmaLine};
    varargout {2} = {DensityTS,WeightedDensityTSSmoothed};
    varargout {3} = {testrfs_lfr_v1v2_freq,testrfs_lfr_v1v2};
elseif nargout == 4 && initialdensitydata == true
    varargout {1} = {PlasmaLine,WeightedPlasmaLine};
    varargout {2} = {DensityTS,WeightedDensityTSSmoothed};
    varargout {3} = {testrfs_lfr_v1v2_freq,testrfs_lfr_v1v2};
    varargout {4} = SubAssistingData;
    
end
end
