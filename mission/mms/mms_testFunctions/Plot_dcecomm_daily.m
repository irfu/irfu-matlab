function Plot_dcecomm_daily(DayOfInterest, bashRun)
% Plot_dcedcv_new is a short function used to plot a full day of all DCE
% and DCV data RAW from our source cdf files. It uses the latest cdf file
% version downloaded from SDC and comine data for an entire day into one
% TSeries. The plot created contains all four sc for quick comparsions.
%
%   Usage: Plot_dcedcv_new('2015-05-20', false);
%   input:  DayOfInterest - string in format 'YYYY-MM-DD'
%           [bashRun] - optional boolean to run in bash mode (save plots to
%           "$DATA_PATH_ROOT/irfu/plots" if it exist).
%
% Created: 2015/06/26
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
% SPDX-License-Identifier: CC-BY-NC-SA-4.0
% See also PLOT_HK10E

narginchk(1,2); % DayOfInterest, [bashRun]
nargoutchk(0,0); % No outputs returned

nowStr = char(datetime("now", "TimeZone","UTC", "Format","uuuu-MM-dd"));

DayOfInterest = cell2mat(regexp(DayOfInterest, '\d{4,4}-\d{2,2}-\d{2,2}', 'match'));
if isempty(DayOfInterest)
  error('Incorrect day string, format should be "YYYY-MM-DD".');
end
if nargin==1, bashRun=false; end % Default to GUI run
dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms/"
outPath = '/data/mms/irfu/plots/edp/RawData/';
if ( ~islogical(bashRun) || ...
    ( bashRun && ~exist(outPath, 'dir') ) )
  error('Incorrect usage. See help.');
end
YLimits.dce = [-15, 15];
YLimits.dcv = [-50, 15];

mms.db_init('local_file_db',dataPathRoot);

% Set default colors as used by Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

% start time is set to midnight of DayOfInterest and stop time end of day.
% NOTE: This quick function does not care about leap second (ie the second
% during 23.59:60.0->23.59.60.999999999Z).
Tint = irf.tint([DayOfInterest,'T00:00:00.000000000Z'], ...
  [DayOfInterest,'T23:59:59.999999999Z']);  %#ok<NASGU>

%% Identify and load all comm files
disp(['Start looking for data from ', DayOfInterest]);
Dce_fast_ts_1=[]; Dce_fast_ts_2=[]; Dce_fast_ts_3=[]; Dce_fast_ts_4=[]; %#ok<NASGU>
Dcv_fast_ts_1=[]; Dcv_fast_ts_2=[]; Dcv_fast_ts_3=[]; Dcv_fast_ts_4=[]; %#ok<NASGU>
Dce_slow_ts_1=[]; Dce_slow_ts_2=[]; Dce_slow_ts_3=[]; Dce_slow_ts_4=[]; %#ok<NASGU>
Dcv_slow_ts_1=[]; Dcv_slow_ts_2=[]; Dcv_slow_ts_3=[]; Dcv_slow_ts_4=[]; %#ok<NASGU>
dce_ts_1=[]; dce_ts_2=[]; dce_ts_3=[]; dce_ts_4=[];
dcv_ts_1=[]; dcv_ts_2=[]; dcv_ts_3=[]; dcv_ts_4=[];
% Create TSeries
c_eval('Dce_fast_ts_?=mms.db_get_ts(''mms?_edp_fast_l1b_dce'', ''mms?_edp_dce_sensor'', Tint);', 1:4);
c_eval('Dcv_fast_ts_?=mms.db_get_ts(''mms?_edp_fast_l1b_dce'', ''mms?_edp_dcv_sensor'', Tint);', 1:4);
c_eval('Dce_slow_ts_?=mms.db_get_ts(''mms?_edp_slow_l1b_dce'', ''mms?_edp_dce_sensor'', Tint);', 1:4);
c_eval('Dcv_slow_ts_?=mms.db_get_ts(''mms?_edp_slow_l1b_dce'', ''mms?_edp_dcv_sensor'', Tint);', 1:4);

emptyFast=false; emptySlow=false;

% Combine (sort data based on time)
for ii=1:4
  c_eval('emptyFast=isempty(Dce_fast_ts_?);', ii);
  c_eval('emptySlow=isempty(Dce_slow_ts_?);', ii);
  if ~emptyFast && ~emptySlow
    % Fast and slow data available (combine them)
    c_eval('dce_ts_?=combine(Dce_fast_ts_?,Dce_slow_ts_?);', ii);
    c_eval('dcv_ts_?=combine(Dcv_fast_ts_?,Dcv_slow_ts_?);', ii);
  elseif ~emptyFast && emptySlow %#ok<UNRCH>
    % Fast data only
    c_eval('dce_ts_?=Dce_fast_ts_?;', ii);
    c_eval('dcv_ts_?=Dcv_fast_ts_?;', ii);
  elseif emptyFast && ~emptySlow
    % Slow data only
    c_eval('dce_ts_?=Dce_slow_ts_?;', ii);
    c_eval('dcv_ts_?=Dcv_slow_ts_?;', ii);
  else
    % Fill with one NaN point to ensure plot function has something to plot.
    c_eval('dce_ts_?=TSeries(EpochTT([DayOfInterest,''T00:00:00.000000000Z'']),NaN);', ii);
    c_eval('dcv_ts_?=TSeries(EpochTT([DayOfInterest,''T00:00:00.000000000Z'']),NaN);', ii);
  end
end

%% Plot All results.
disp('Plotting all data from all S/C.');

fig.dce = figure; h.dce = irf_plot({dce_ts_1, dce_ts_2, dce_ts_3, dce_ts_4});
fig.dcv = figure; h.dcv = irf_plot({dcv_ts_1, dcv_ts_2, dcv_ts_3, dcv_ts_4});

title(h.dce(1),['Plot created: ', nowStr, '. DCE raw for all four s/c.']);
legend(h.dce(1), 'DCE 12', 'DCE 34', 'DCE 56');

title(h.dcv(1),['Plot created: ', nowStr, '. DCV raw for all four s/c.']);
legend(h.dcv(1), 'DCV 1', 'DCV 3', 'DCV 5');

% Improve readability
% Zoom in to aviod Sweep spikes (typically probe potential should be less
% than this absolute value in normal operations).
rescaleY(h.dce,YLimits.dce); rescaleY(h.dcv,YLimits.dcv);
for ii=1:4
  % Split YLabel to use two lines with "units" on the second line.
  ylabel(h.dce(ii), sprintf('MMS%i\n[mV/m]', ii));
  ylabel(h.dcv(ii), sprintf('MMS%i\n[V]', ii));
end

if(bashRun)
  disp('Saving entire day');
  iData = {'dce', 'dcv'};
  for ii=1:length(iData)
    outFolder=[outPath, iData{ii}, filesep, DayOfInterest(1:4)];
    % Create folder (if new year)
    if ~exist(outFolder,'dir'), mkdir(outFolder); end
    filePrefix = [DayOfInterest, '_', iData{ii}, '_'];
    % Entire day:
    % savefig(fig.(iData{ii}), [outPath, DayOfInterest, '_', iData{ii}, '_000000_235959'], 'compact'); %R2013b or later
    print(fig.(iData{ii}), '-dpng', [outFolder, filesep, filePrefix, '000000_235959']);

    % Interval 00:01:00 to 00:11:00 UTC
    tint = irf.tint([DayOfInterest,'T00:01:00Z'], [DayOfInterest,'T00:11:00Z']);
    irf_zoom(h.(iData{ii}), 'x', tint);
    irf_zoom(h.(iData{ii}), 'y');
    rescaleY(h.(iData{ii}), YLimits.(iData{ii}));
    title(h.(iData{ii})(1),{['Plot created: ', nowStr, '. ', upper(iData{ii}), ' raw for all four s/c.'],...
      'Zoom in: 00:01:00 to 00:11:00 UTC'});
    disp('Saving 00:01:00 to 00:11:00');
    print(fig.(iData{ii}), '-dpng', [outFolder, filesep, filePrefix, '000100_001100']);

    % Interval 06:00:00 to 06:10:00 UTC
    tint = irf.tint([DayOfInterest,'T06:00:00Z'], [DayOfInterest,'T06:10:00Z']);
    irf_zoom(h.(iData{ii}), 'x', tint);
    irf_zoom(h.(iData{ii}),'y');
    rescaleY(h.(iData{ii}), YLimits.(iData{ii}));
    title(h.(iData{ii})(1),{['Plot created: ',nowStr,'. ', upper(iData{ii}), ' raw for all four s/c.'],...
      'Zoom in: 06:00:00 to 06:10:00 UTC.'});
    disp('Saving 06:00:00 to 06:10:00');
    print(fig.(iData{ii}), '-dpng', [outFolder, filesep, filePrefix, '060000_061000']);

    % Interval 12:00:00 to 12:01:00 UTC (ie a shorter interval as well).
    tint = irf.tint([DayOfInterest,'T12:00:00Z'], [DayOfInterest,'T12:01:00Z']);
    irf_zoom(h.(iData{ii}), 'x', tint);
    irf_zoom(h.(iData{ii}),'y');
    rescaleY(h.(iData{ii}), YLimits.(iData{ii}));
    title(h.(iData{ii})(1),{['Plot created: ',nowStr,'. ', upper(iData{ii}), ' raw for all four s/c.'],...
      'Zoom in: 12:00:00 to 12:01:00 UTC.'});
    disp('Saving 12:00:00 to 12:01:00');
    print(fig.(iData{ii}), '-dpng', [outFolder, filesep, filePrefix, '120000_120100']);

    % Interval 18:00:00 to 18:10:00 UTC
    tint = irf.tint([DayOfInterest,'T18:00:00Z'], [DayOfInterest,'T18:10:00Z']);
    irf_zoom(h.(iData{ii}), 'x', tint);
    irf_zoom(h.(iData{ii}), 'y');
    rescaleY(h.(iData{ii}), YLimits.(iData{ii}));
    title(h.(iData{ii})(1),{['Plot created: ',nowStr,'. ', upper(iData{ii}), ' raw for all four s/c.'],...
      'Zoom in: 18:00:00 to 18:10:00 UTC.'});
    disp('Saving 18:00:00 to 18:10:00');
    print(fig.(iData{ii}), '-dpng', [outFolder, filesep, filePrefix, '180000_181000']);
  end
end

  function rescaleY(h,MinMax)
    % Rescale h so its max/min is +/-15 or less.
    narginchk(2,2);
    for id=1:4
      % Get currently used limits
      ylimits = ylim(h(id));
      % If outside region -15 to +15, rescale.
      if(ylimits(1)<-MinMax(1) || ylimits(2)>MinMax(2))
        ylimits(1) = max(MinMax(1), ylimits(1));
        ylimits(2) = min(MinMax(2), max(ylimits(2),MinMax(1)+1));
        ylim(h(id), ylimits);
      end
    end
  end

end
