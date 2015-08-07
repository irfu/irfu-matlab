function Plot_mms_results(Day, varargin)
% Plot_MMS_RESULTS function for plotting data from MMS in order to locate
% any issues. Will read all data avaiable for each day.
%
%  Inputs:
%    DayOfInterest - day string in format 'YYYY-MM-DD'.
%    [bashRun]     - optional boolean to run in bash mode (save plots to
%           "$DATA_PATH_ROOT/irfu/plots" if it exist).
%  Usage:
%   Plot_mms_results('2015-05-12','bashRun',true,'plotUsc',false);
%  to automatically generate plots for all but Usc.
%
% Spectra (of DSL E-x & E-y for all MMS s/c) will cover one figure (if it's
% to be plotted).
% Usc (of all MMS s/c) will be plotted on a second figure (if it's to be
% plotted).
% DSL E-x & E-y (for all MMS s/c) will also be plotted on this second
% figure (if it's to be plotted).
%
% Created: 2015/08/04
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
%

% Process input arguments
p = inputParser;
% Default values
default.BashRun = true; % Assume automated bashRun
default.PlotUsc = true; % Assume Usc should be plotted.
default.PlotDSL = true; % Assume DSL E-x & E-y should be plotted.
default.PlotSpect = true; % Assume DSL E-x & E-y spectra should be plotted.
validDay = @(x) assert(ischar(x) && ...
  numel(sscanf(x,'%4d-%2d-%2d%s'))==3 && ...
  length(x)==10, ...
  'Format should be "YYYY-MM-DD".'); % Day must be correct format and char.
% Input arguments, processed in this order of not given as explicit
% arguments or specified in a struct.
addRequired(p, 'Day', validDay);
addOptional(p, 'bashRun', default.BashRun, @islogical);
addOptional(p, 'plotUsc', default.PlotUsc, @islogical);
addOptional(p, 'plotDSL', default.PlotDSL, @islogical);
addOptional(p, 'plotSpect', default.PlotSpect, @islogical);

parse(p, Day, varargin{:});



nowStr = irf_time(now,'datenum>utc_yyyy-mm-dd');

dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms/"
mms.db_init('local_file_db',dataPathRoot);
MMS_CONST = mms_constants();
% Set default colors as used by Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

if(p.Results.bashRun)
  outPath = [dataPathRoot,'irfu',filesep,'plots',filesep];
  if(~exist(outPath,'dir')), error('outpath does not exist'); end
end

% Create time interval covering the entire day. Note, this will skip the
% leap second (interval such as 2015-07-01T23:59:60.000000000Z to
% 2015-07-01T23:59.60.999999999Z).
tint = irf.tint([p.Results.Day,'T00:00:00.000000000Z'],...
  [p.Results.Day,'T23:59:59.999999999Z']);

SCid = {'mms1', 'mms2', 'mms3', 'mms4'};
tmMode = {'comm', 'fast', 'slow', 'brst'};


%% Identify and load all data files.
disp(['Start looking data from ', p.Results.Day]);

% File suffix
fileSuf = {'l2_scpot','ql_dce2d'};
% Data suffix
dataSuf = {'scpot','dce_xyz_dsl'};

% Identify and load all "comm", "slow", "fast" and "brst" mode files.
if(p.Results.plotUsc)
  uscTs = getDATA(1); % ALL Usc data for the day as struct of a TSeries
end


if(p.Results.plotDSL)
  dslTs = getDATA(2); % ALL DSL data for the day as struct of a TSeries
end



if(p.Results.plotSpect)
  % Get only comm data (Almost exclusively 128 Hz data)
  tmMode = {'comm'};
  dslCommTs = getDATA(2); % DSL data of only comm ("128Hz") as a struct of TSeries

  %% Compute power spectra of DSL data
  % Comm 128 Hz
  for kk=1:length(SCid)
    disp('Computing spectra for ',SCid{kk});
    dslCommSpect.(SCid{kk}) = irf_powerfft(dslCommTs.(SCid{kk}), ...
      8192, ...
      MMS_CONST.Samplerate.comm_128);
    if(isempty(dslCommSpect.(SCid{kk})))
      % FIXME better..
      warning(['Empty spectra from DSL on ', SCid{kk}]);
    end
  end

  % Or only DSL E X component
  %dslCommSpect.(SCid{kk}) = irf_powerfft(dslCommTs.(SCid{kk}).x, ...
  %  8192, ...
  %  MMS_CONST.Samplerate.comm_128);


  % Get fast (32Hz) or slow (8Hz) mode data, starting around 2015-07-16
  %tmMode = {'fast'};
  %dslFastTs = getDATA(2); % dsl data of only fast as a TSeries

  %tmMode = {'slow'};
  %dslSlowTs = getDATA(2); % dsl data of only slow as a TSeries

  %dslFastSpect = irf_powerfft(dslFastTs.data, ...
  %  512, ...
  %  MMS_CONST.Samplerate.fast);
  %dslFastSpect = irf_powerfft(dslSlowTs.data, ...
  %  512, ...
  %  MMS_CONST.Samplerate.slow);

  % FIXME:
  % Combine the various power spectra in some nice function so plots are
  % continuous...

  % Restore tmModes
  tmMode = {'comm','fast','slow','brst'};

end





%% Plot data

disp('Plotting all data from all S/C.');

if(p.Results.plotUsc || p.Results.plotDSL)
  % Create a new figure, with 3 subplots. First will be DSL E-x, then DSL
  % E-y and third Usc. (Regardless if DSL is to be plotted? Keep position
  % in figure to ensure easy comparison between pictures?)
  h = irf_plot(3, 'newfigure'); % FIXME, see below... (how to plot to specific subplots.)

  if(p.Results.plotUsc)
    % PLOT USC
    disp('Plotting Usc');
    % toPlot = cell(1,length(SCid));
    % for i=1:length(SCid)
    %  toPlot{1,i} = uscTs.(SCid{i});
    % end
    % h = irf_plot(toPlot,'comp'); % FIXME, plot only to subplot 1 (of h).
    % title(h(1),['Plot created: ',nowStr,'. Usc from SDC files for all four s/c.']);
    % legend(h(1), 'MMS1', 'MMS2', 'MMS3', 'MMS4');
    % ylabel(h(1),{'Spacecraft potential','[V]'});
  end

  if(p.Results.plotDSL)
    % PLOT DSL E-X, E-Y
    disp('Plotting DSL');
    toPlot = cell(1,length(SCid));
    for i=1:length(SCid)
      toPlot{1,i} = dslTs.(SCid{i});
    end
    h = irf_plot(toPlot,'comp'); % FIXME, plot only DSL E-x & E-y in the
    % subplot 2 and 3 (of h) created above just after beginning of plotUsc or
    % plotDSL
    %title(h(1),['Plot created: ',nowStr,'. DSL from SDC files for all four s/c.']);
    legend(h(1), 'MMS1', 'MMS2', 'MMS3', 'MMS4');
    ylabel(h(1),{'DSL E-x','[mV/m]'});
    ylabel(h(2),{'DSL E-y','[mV/m]'});
    %ylabel(h(3),{'BCS E-z','[mV/m]'});
  end

  if(p.Results.bashRun)
    % Save plots
  end
end % End of plot figure for Usc and/or DSL E-x & E-y



if(p.Results.plotSpect)
  % PLOT spectogram of DSL
  disp('Plotting DSL spectra');

  yTicks = [0.05, 0.10, 0.15, 0.2, 0.4, 0.8, 1.0, 10];
  yLimits = [0, 12.5];
  h = irf_plot(8, 'newfigure');

  h = irf_spectrogram(dslCommSpect.mms1); % FIXME plot to specific subplot (of h).
  h(1).YScale = 'log';
  h(1).YTick = yTicks;
  h(1).YLim = yLimits;
  ylabel(h(1),{'MMS1 DSL E-x','freq [Hz]'});
  caxis(h(1),[-4 1]);

  % if multiple components:
  h(2).YScale = 'log';
  h(2).YTick = yTicks;
  h(2).YLim = yLimits;
  ylabel(h(2),{'MMS1 DSL E-y','freq [Hz]'});
  caxis(h(2),[-4 1]);

  % FIXME, should perhaps not be plotted..
  h(3).YScale = 'log';
  h(3).YTick = yTicks;
  h(3).YLim = yLimits;
  ylabel(h(3),{'MMS1 BCS E-z','freq [Hz]'});
  caxis(h(3),[-4 1]);

  % fig = figure;
  % h = irf_spectrogram(struct('f',dslCommSpect.mms1.f,'t',...
  %   dslCommSpect.mms1.t,'p',dslCommSpect.mms1.p(2)));
  % h(1).YScale = 'log';
  % h(1).YTick = [.25, 0.5, 1, 10]; % Limits from Cluster, may need updating.
  % h(1).YLim = [0, 12.5];
  % ylabel(h(1),{'MMS1 DSL E-y power spectra','f [Hz]'});

  if(p.Results.bashRun)
    % Save plots
  end
end % End of plot figure with DSL spectra





%% Help function to locate and load MMS data into TSeries
function outTs = getDATA(dataID)
  if(dataID>size(dataSuf)), error('Incorrect usage'); end;
  for ii=1:length(SCid)
    disp(['Loading and converting data from ',SCid{ii},' for ',dataSuf{dataID}]);

    TsCombined = struct(SCid{ii},struct('time',[],'data',[]));

    for jj=1:length(tmMode)
      TsTmp.(SCid{ii}).(tmMode{jj}) = mms.db_get_ts([SCid{ii},'_edp_',tmMode{jj},'_',fileSuf{dataID}],...
        [SCid{ii},'_edp_',dataSuf{dataID}], tint);
      if(isempty(TsTmp.(SCid{ii}).(tmMode{jj})))
        warning(['No ',SCid{ii},' ', tmMode{jj}, ' ', dataSuf{dataID},...
          ' file was found for this day.']);
        continue
      end
      % Combine the data and time series
      TsCombined.(SCid{ii}).time = [TsCombined.(SCid{ii}).time; TsTmp.(SCid{ii}).(tmMode{jj}).time.ttns];
      TsCombined.(SCid{ii}).data = [TsCombined.(SCid{ii}).data; TsTmp.(SCid{ii}).(tmMode{jj}).data];
    end

    % Sort by time
    [~, idxSort] = sort(TsCombined.(SCid{ii}).time);
    TsCombined.(SCid{ii}).time = TsCombined.(SCid{ii}).time(idxSort);
    TsCombined.(SCid{ii}).data = TsCombined.(SCid{ii}).data(idxSort,:);
    % Remove duplicates
    [~,idxUniq] = unique(TsCombined.(SCid{ii}).time);
    TsCombined.(SCid{ii}).time = TsCombined.(SCid{ii}).time(idxUniq);
    TsCombined.(SCid{ii}).data = TsCombined.(SCid{ii}).data(idxUniq,:);

    % Return a struct with fields of time series for each sc.
    if(~isempty(TsCombined.(SCid{ii}).time))
      outTs.(SCid{ii}) = TSeries(EpochTT(TsCombined.(SCid{ii}).time),TsCombined.(SCid{ii}).data);
    else
      % Fill with one NaN data point if no data was read at all for this sc.
      outTs.(SCid{ii}) = TSeries(EpochTT(tint.start), NaN(1,3), 'vec_xyz');
    end
    outTs.(SCid{ii}).name = upper(SCid{ii});
  end
end

end