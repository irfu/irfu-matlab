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
% Or use a struct:
%  inArg = struct('bashRun', false, 'plotUsc', false, 'plotDSL', false, ...
%            'plotSpect', true);
%  Plot_mms_results('2015-05-12', inArg)
%
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
% SPDX-License-Identifier: CC-BY-NC-SA-4.0
tic;
p = verify_input(Day, varargin{:});

yyyy = str2double(Day(1:4)); % Store numerical value
mm = str2double(Day(6:7)); dd = str2double(Day(9:10));

nowStr = irf_time(now,'datenum>utc_yyyy-mm-dd');

dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms"
mms.db_init('local_file_db',dataPathRoot);
MMS_CONST = mms_constants();
% Set default colors as used by Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

if(p.Results.bashRun)
  outPath = [dataPathRoot,'irfu',filesep,'plots',filesep,'results',filesep];
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
  if(tint.start<EpochTT('2015-07-28T00:00:00')) % Last day when any Comm files existed
    % Assume it is only Comm data
    %   % Get only comm data (Almost exclusively 128 Hz data)
    tmMode = {'comm'};
    dslCommTs = getDATA(2); % DSL data of only comm ("128Hz") as a struct of TSeries
    %
    %   %% Compute power spectra of DSL data
    %   % Comm 128 Hz
    for kk=1:length(SCid)
      disp(['Computing spectra for DSL xyz for ',SCid{kk}]);
      dslComm128Spect.(SCid{kk}) = irf_powerfft(dslCommTs.(SCid{kk}), ...
        8192, ...
        MMS_CONST.Samplerate.comm_128);
      if(isempty(dslComm128Spect.(SCid{kk})))
        % FIXME better..
        warning(['Empty spectra from DSL on ', SCid{kk}]);
      end
      dslCommTs.(SCid{kk})=[]; % Empty to save memory.
    end
    
  end
  
  % Get fast (32Hz) or slow (8Hz) mode data, starting around 2015-07-15
  tmMode = {'fast'};
  dslFastTs = getDATA(2); % dsl data of only fast as a TSeries
  
  tmMode = {'slow'};
  dslSlowTs = getDATA(2); % dsl data of only slow as a TSeries
  
  for kk=1:length(SCid)
    disp(['Computing fast spectra for DSL xyz for ',SCid{kk}]);
    dslFastSpect.(SCid{kk}) = irf_powerfft(dslFastTs.(SCid{kk}), ...
      512, ...
      MMS_CONST.Samplerate.fast);
    if(isempty(dslFastSpect.(SCid{kk})))
      % FIXME better..
      warning(['Empty spectra from DSL on ', SCid{kk}]);
    end
    dslFastTs.(SCid{kk})=[]; % Empty to save memory.
    
    disp(['Computing slow spectra for DSL xyz for ',SCid{kk}]);
    slowRate = MMS_CONST.Samplerate.slow{1};
    dtSlow = round(10^9/median(double(diff(dslSlowTs.(SCid{kk}).time.ttns))));
    for iSmplrate=1:length(MMS_CONST.Samplerate.slow)
      if(MMS_CONST.Samplerate.slow{iSmplrate}*0.9 <= dtSlow && dtSlow <= MMS_CONST.Samplerate.slow{iSmplrate}*1.1)
        slowRate = MMS_CONST.Samplerate.slow{iSmplrate};
      end
    end
    dslSlowSpect.(SCid{kk}) = irf_powerfft(dslSlowTs.(SCid{kk}), ...
      128, ...
      slowRate);
    if(isempty(dslSlowSpect.(SCid{kk})))
      % FIXME better..
      warning(['Empty spectra from DSL on ', SCid{kk}]);
    end
    dslSlowTs.(SCid{kk})=[]; % Empty to save memory.
    
    % Combine the various power spectra in some nice function so plots are
    % continuous...
    % THIS WORKS but will slightly adjust times.
    comb.(SCid{kk}) = struct('f',0:0.0625:15.9375);
    newTime = tint.start:15.9688:tint.stop;
    newTime = newTime.epochUnix;
    comb.(SCid{kk}).t = newTime;
    len_time = size(newTime,1);
    % len_freq = size(dslFastSpect.(SCid{kk}).f,1) - size(dslSlowSpect.(SCid{kk}).f,1);
    len_freq = 192; % = 256 - 64, size(Fast.f) - slow(Slow.f)
    for ll=1:3 % 3 dimensions (DSL x, y, z)
      comb.(SCid{kk}).p{1,ll} = NaN(len_time, 256);
      if(~isempty(dslFastSpect.(SCid{kk})))
        if(length(dslFastSpect.(SCid{kk}).t)>1)
          dslFastNew.p{1,ll} = (interp2((dslFastSpect.(SCid{kk}).t)',...
            dslFastSpect.(SCid{kk}).f, (dslFastSpect.(SCid{kk}).p{1,ll})', ...
            newTime', dslFastSpect.(SCid{kk}).f))';
          indNotNaN = ~all(isnan(dslFastNew.p{1,ll}),2);
          comb.(SCid{kk}).p{1,ll}(indNotNaN,:) = dslFastNew.p{1,ll}(indNotNaN,:);
        end
      end
      if(~isempty(dslSlowSpect.(SCid{kk})))
        if(length(dslSlowSpect.(SCid{kk}).t)>1)
          dslSlowNew.p{1,ll} = (interp2((dslSlowSpect.(SCid{kk}).t)',...
            dslSlowSpect.(SCid{kk}).f, (dslSlowSpect.(SCid{kk}).p{1,ll})', ...
            newTime', dslSlowSpect.(SCid{kk}).f))';
          % Fill high frequencies (that only exist in Fast) with NaN in the slow spectra
          dslSlowNew.p{1,ll} = [dslSlowNew.p{1,ll} NaN(len_time, len_freq)];
          indNotNaN = ~all(isnan(dslSlowNew.p{1,ll}),2);
          comb.(SCid{kk}).p{1,ll}(indNotNaN,:) = dslSlowNew.p{1,ll}(indNotNaN,:);
        end
      end
      if(tint.start<EpochTT('2015-07-28T00:00:00') && ~isempty(dslComm128Spect.(SCid{kk})))
        dslCommNew.p{1,ll} = (interp2((dslComm128Spect.(SCid{kk}).t)',...
          dslComm128Spect.(SCid{kk}).f, (dslComm128Spect.(SCid{kk}).p{1,ll})',...
          newTime', (0:0.0625:15.9375)'))';
        indNotNaN = ~all(isnan(dslCommNew.p{1,ll}),2);
        comb.(SCid{kk}).p{1,ll}(indNotNaN,:) = dslCommNew.p{1,ll}(indNotNaN,:);
      end
    end
    
    dslFastSpect.(SCid{kk}) = []; % Save some memory
    dslComm128Spect.(SCid{kk}) = [];
    dslSlowSpect.(SCid{kk}) = [];
    
  end
  dslCommSpect = comb;
  clear comb newTime dslSlowTs dslFastTs dslSlowNew dslFastNew;
  clear dslSlowSpect dslFastSpect dslComm128Spect;
  %  end
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
    toPlot = cell(1,length(SCid));
    for i=1:length(SCid)
      toPlot{1,i} = uscTs.(SCid{i});
    end
    irf_plot(h(1),toPlot,'comp');
    title(h(1),['Plot created: ',nowStr,'. Usc and DSL X&Y.']);
    legend(h(1), 'MMS1', 'MMS2', 'MMS3', 'MMS4');
    ylabel(h(1),{'Spacecraft potential','[V]'});
    clear toPlot uscTs; % save some memory...
  end
  
  if(p.Results.plotDSL)
    % PLOT DSL E-X, E-Y
    disp('Plotting DSL');
    toPlot = cell(2,length(SCid));
    for i=1:length(SCid)
      toPlot{1,i} = dslTs.(SCid{i}).x;
      toPlot{2,i} = dslTs.(SCid{i}).y;
    end
    irf_plot(h(2),toPlot(1,:),'comp');
    irf_plot(h(3),toPlot(2,:),'comp');
    % Zoom in
    irf_zoom(h(2),'y',[-10 10])
    irf_zoom(h(3),'y',[-10 10])
    legend(h(1), 'MMS1', 'MMS2', 'MMS3', 'MMS4');
    ylabel(h(2),{'DSL E-x','[mV/m]'});
    ylabel(h(3),{'DSL E-y','[mV/m]'});
    clear toPlot dslTs; % save some memory..
  end
  
  if(p.Results.bashRun)
    % Save plots in 3 hour segments, ie 00:00:00 to 03:00:00 and so on
    for iZoom=0:3:21
      tintZoom = {[yyyy, mm, dd, iZoom, 00, 00], ...
        [yyyy, mm, dd, iZoom+3, 00, 00]};
      irf_zoom(h,'x',tintZoom);
      % Save segment
      print(gcf, '-dpng', sprintf('%s%sT%02d0000_%02d0000_Usc_and_DSLxy',...
        outPath, Day, iZoom, iZoom+3) );
      toc;
    end
  end
end % End of plot figure for Usc and/or DSL E-x & E-y



if(p.Results.plotSpect)
  clear dslCommTs
  % PLOT spectogram of DSL
  disp('Plotting DSL spectra');
  
  yTicks = [0.1, 1.0, 10];
  yTickLabel = {0.1; 1.0; 10};
  yLimits = [0, 12.5];
  h = irf_plot(length(SCid), 'newfigure');
  set(gcf,'position',[7   159   790   916]);
  figure_start_epoch(tint.start.epochUnix);
  
  for kk=1:length(SCid)
    % DSL E-X
    hold(h(kk),'on')
    % if both X & Y spectra was computed and is to be added into one
    % subplot, then use:
    irf_spectrogram(h(kk), struct('f', dslCommSpect.(SCid{kk}).f, ...
      't', dslCommSpect.(SCid{kk}).t, 'p', dslCommSpect.(SCid{kk}).p{1} + ...
      dslCommSpect.(SCid{kk}).p{2}));
    % Create spectrogram of only DSL X:
    %irf_spectrogram(h(kk), dslCommSpect.(SCid{kk}));
    set(h(kk),'YTick',yTicks,'YScale','log', 'YTickLabel',yTickLabel);
    h(kk).YLim = yLimits;
    ylabel(h(kk),{upper(SCid{kk}),'freq [Hz]'});
    caxis(h(kk),[-4 1]);
    hold(h(kk),'off')
    if(kk==1), title(h(kk),['Plot created: ',nowStr, ...
        '. Combined pds(Ex_{dsl}) + pds(Ey_{dsl}).']); end
    dslCommSpect.(SCid{kk}) = []; % Clear to save memory
  end
  
  clear dslCommSpect
  
  if(p.Results.bashRun)
    % Save plots in 3 hour segments, ie 00:00:00 to 03:00:00 and so on
    for iZoom=0:3:21
      tintZoom = {[yyyy, mm, dd, iZoom, 00, 00], ...
        [yyyy, mm, dd, iZoom+3, 00, 00]};
      irf_zoom(h,'x',tintZoom);
      % Save segment
      print(gcf, '-dpng', sprintf('%s%sT%02d0000_%02d0000_DSLxy_spectra',...
        outPath, Day, iZoom, iZoom+3) );
      toc;
    end
  end
end % End of plot figure with DSL spectra





%% Help function to locate and load MMS data into TSeries
  function outTs = getDATA(dataID)
    if(dataID>size(dataSuf)), error('Incorrect usage'); end
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
        if(size(TsCombined.(SCid{ii}).data,2)==3)
          outTs.(SCid{ii}) = TSeries(EpochTT(TsCombined.(SCid{ii}).time),TsCombined.(SCid{ii}).data,'vec_xyz');
        else
          outTs.(SCid{ii}) = TSeries(EpochTT(TsCombined.(SCid{ii}).time),TsCombined.(SCid{ii}).data);
        end
      else
        % Fill with one NaN data point if no data was read at all for this sc.
        outTs.(SCid{ii}) = TSeries(EpochTT(tint.start), NaN(1,3), 'vec_xyz');
      end
      outTs.(SCid{ii}).name = upper(SCid{ii});
      
      clear TsCombined TsTmp % Clear to save memory
    end
  end

  function p = verify_input(Day, tmpVarargin)
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
    parse(p, Day, tmpVarargin);
  end

  function t_start_epoch = figure_start_epoch(st)
    ud = get(gcf,'userdata');
    if isfield(ud,'t_start_epoch')
      t_start_epoch = ud.t_start_epoch;
    else
      ud.t_start_epoch = st;
      set(gcf,'userdata',ud);
    end
  end

end
