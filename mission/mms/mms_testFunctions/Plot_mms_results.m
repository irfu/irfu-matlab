function Plot_mms_results(DayOfInterest, bashRun)
% Plot_MMS_RESULTS function for plotting data from MMS in order to locate
% any issues. Will read all data avaiable for each day.
%
%  Inputs:
%    DayOfInterest - day string in format 'YYYY-MM-DD'.
%    [bashRun]     - optional boolean to run in bash mode (save plots to
%           "$DATA_PATH_ROOT/irfu/plots" if it exist).
%  Usage:
%   Plot_mms_results('2015-05-12',true)
%
% Created: 2015/08/04
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
%

narginchk(1,2);
if(nargin==1), bashRun = false; end

nowStr = irf_time(now,'datenum>utc_yyyy-mm-dd');
if((numel(sscanf(DayOfInterest,'%4d-%2d-%2d%s'))~=3)||(length(DayOfInterest)~=10))
  error('Incorrect day string, format should be "YYYY-MM-DD".');
end
dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms/"
mms.db_init('local_file_db',dataPathRoot);
MMS_CONST = mms_constants();
% Set default colors as used by Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

if(bashRun)
  outPath = [dataPathRoot,'irfu',filesep,'plots',filesep];
  if(~exist(outPath,'dir')), error('outpath does not exist'); end
end

% Create time interval covering the entire day. Note, this will skip the
% leap second (interval such as 2015-07-01T23:59:60.000000000Z to
% 2015-07-01T23:59.60.999999999Z).
tint = irf.tint([DayOfInterest,'T00:00:00.000000000Z'],...
  [DayOfInterest,'T23:59:59.999999999Z']);

SCid = {'mms1', 'mms2', 'mms3', 'mms4'};
tmMode = {'comm', 'fast', 'slow', 'brst'};

%% Identify and load all data files.
disp(['Start looking data from ', DayOfInterest]);

% File suffix
fileSuf = {'l2_scpot','ql_dce2d'};
% Data suffix
dataSuf = {'scpot','dce_xyz_dsl'};
% Identify and load all "comm", "slow", "fast" and "brst" mode files.
uscTs = getDATA(1);
%dslTs = getDATA(2); % ALL dsl data as a TSeries

% Get only comm data (Almost exclusively 128 Hz data)
tmMode = {'comm'};
dslCommTs = getDATA(2); % dsl data of only comm "128Hz" as a TSeries

% Get fast (32Hz) or slow (8Hz) mode data, starting around 2015-07-16
%tmMode = {'fast'};
%dslFastTs = getDATA(2); % dsl data of only fast as a TSeries
%tmMode = {'slow'};
%dslSlowTs = getDATA(2); % dsl data of only slow as a TSeries


% Reset tmModes
tmMode = {'comm','fast','slow','brst'};

%% Compute power spectra of DSL data
% Comm 128 Hz
dslCommSpect.mms1 = irf_powerfft(dslCommTs.mms1, ...
  8192, ...
  MMS_CONST.Samplerate.comm_128);


%dslFastSpect = irf_powerfft(dslFastTs.data, ...
%  512, ...
%  MMS_CONST.Samplerate.fast);
%dslFastSpect = irf_powerfft(dslSlowTs.data, ...
%  512, ...
%  MMS_CONST.Samplerate.slow);

% FIXME:
% Combine the various power spectra in some nice function so plots are
% continuous...


%% Plot data

disp('Plotting all data from all S/C.');

% PLOT USC
toPlot = cell(1,length(SCid));
for i=1:length(SCid)
  toPlot{1,i} = uscTs.(SCid{i});
end
fig = figure;
h = irf_plot(toPlot,'comp');
title(h(1),['Plot created: ',nowStr,'. Usc from SDC files for all four s/c.']);
legend(h(1), 'MMS1', 'MMS2', 'MMS3', 'MMS4');
ylabel(h(1),{'Spacecraft potential','[V]'});


% Plot spectogram of DSL
fig = figure;
h = irf_spectrogram(dslCommSpect.mms1);
h.YScale = 'log';
h.YTick = [.25, 0.5, 1, 10]; % Limits from Cluster, may need updating.
h.YLim = [0, 12.5];
ylabel(h,{'MMS1','f [Hz]'});



% Help function to locate and load MMS data into TSeries
function outTs = getDATA(dataID)
  if(dataID>size(dataSuf)), error('Incorrect usage'); end;
  for ii=1:length(SCid)
    disp(['Loading and converting data from ',SCid{ii},' for ',dataSuf{dataID}]);

    TsCombined = struct(SCid{ii},struct('time',[],'data',[]));

    for jj=1:length(tmMode)
      TsTmp.(SCid{ii}).(tmMode{jj}) = mms.db_get_ts([SCid{ii},'_edp_',tmMode{jj},'_',fileSuf{dataID}],...
        [SCid{ii},'_edp_',dataSuf{dataID}], tint);
      if(isempty(TsTmp.(SCid{ii}).(tmMode{jj})))
        warning(['No ',SCid{ii},' ', tmMode{jj}, dataSuf{dataID},' file was found for this day.']);
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
      % Fill with one empty data point if no data was read at all for this sc.
      outTs.(SCid{ii}) = TSeries(EpochTT(tint.start),[0, 0, 0]);
    end
    outTs.(SCid{ii}).name = upper(SCid{ii});
  end
end

end