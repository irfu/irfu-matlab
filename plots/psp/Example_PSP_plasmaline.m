% Example plasma line tracking plot & comparison with SPI

start = [2020 06 09 03 00 00] ;
stop =  [2020 06 09 04 00 00];
fullday = [2020 06 09] ;
SpectrumData = psp_load([],'pl_rfs_lfr', [start(1) start(2) start(3)], [stop(1) stop(2) stop(3)]); %

%% Some examples

% Let's give no initial density data, but a bin, let's say 30 to start and
% a +- 3 bins to search for.
psp_freqtracker(start,stop,'initialdensitydata',false,'initialfreq',30,'binrange',3,'generateplot',true);

% Quite a few spikes, so maybe initial data and leaving the bins at default
% (+-2) will give a better plot
[PlasmaLine,DensityTS,Pspdata,InitialDensityTS] = psp_freqtracker(start,stop,'initialdensitydata',true,'generateplot',true,'densitycomparisonplot',true);

% One can also track a full day but this will generally provide a weaker
% characterization.
%[PlasmaLine,~,Pspdata,~] = psp_freqtracker(fullday,'generateplot',false);

