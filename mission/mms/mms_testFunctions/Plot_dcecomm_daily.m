function [h, dceTs] = Plot_dcecomm_daily(DayOfInterest, bashRun)
% Plot_dce128_daily is a short function used to plot a full day of all DCE
% data RAW from our source cdf files. It uses the latest cdf file version
% downloaded from SDC and comine data for an entire day into one TSeries.
% The plot created contains all four sc for quick and easy comparison.
%
%   Usage: [[h], [dceTs]] = PLOT_DCECOMM_DAILY('2015-05-20');
%   input:  DayOfInterest - string in format 'YYYY-MM-DD'
%           [bashRun] - optional boolean to run in bash mode (save plots to
%           "$DATA_PATH_ROOT/irfu/plots" if it exist).
%   output: [h]     - optional axes handle for irf_plot
%           [dceTs] - optional struct of TSeries dce data.
%
% Created: 2015/06/26
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
%
% See also PLOT_HK10E

narginchk(1,2); % DayOfInterest, [bashRun]
nargoutchk(0,2); % [h], [dceTs] 

nowStr = irf_time(now,'datenum>utc_yyyy-mm-dd');
if((numel(sscanf(DayOfInterest,'%4d-%2d-%2d%s'))~=3)||(length(DayOfInterest)~=10))
  error('Incorrect day string, format should be "YYYY-MM-DD".');
end
if(nargin==1), bashRun=false; end; % Default to GUI run
dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms/"
if( ~islogical(bashRun) || ...
    ( bashRun && ~exist([dataPathRoot,'irfu',filesep,'plots'],'dir') ) )
  error('Incorrect usage. See help.');
end
mms.db_init('local_file_db',dataPathRoot)
% Set default colors as used by Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

% start time is set to midnight of DayOfInterest and stop time end of day.
% NOTE: This quick function does not care about leap second (ie the second
% during 23.59:60.0->23.59.60.999999999Z).
tint = irf.tint([DayOfInterest,'T00:00:00.000000000Z'],[DayOfInterest,'T23:59:59.999999999Z']);

%% Identify and load all DCE_128 files.
SCid = {'mms1', 'mms2', 'mms3', 'mms4'};

for ii=1:length(SCid)
  disp(['Loading and converting data from ',SCid{ii}]);
  % Combine DCE to one TSeries per day and sc
  dceTs.(SCid{ii}) = mms.db_get_ts([SCid{ii},'_edp_comm_l1b_dcecomm'],...
    [SCid{ii},'_edp_dce_sensor'], tint);
  % Note: DC V cannot easily be converted to a TSeries as it sometimes is
  % listed in the dce file and some other times contain all six probes...
  % dcv_ts=mms.db_get_ts('mms1_edp_comm_l1b_dcv128','mms1_edp_dcv_sensor',tint);
  if(isempty(dceTs.(SCid{ii})))
    warning(['No ',SCid{ii},' DCECOMM file found for this day.']);
  end
end

%% Plot All results.
% Note this will currently fail if not at least one DCECOMM files was
% identified for each s/c for the day of interest.
disp('Plotting all data from all S/C.');
toPlot = cell(1,length(SCid));
for ii=1:length(SCid)
  toPlot{1,ii} = dceTs.(SCid{ii});
end
fig = figure;
h = irf_plot(toPlot);
title(h(1),['Plot created: ',nowStr,'. DCE raw from our source cdf files for all four s/c.']);
legend(h(1), 'DCE 12', 'DCE 34', 'DCE 56');

% Improve readability
% Zoom in to aviod Sweep spikes (typically probe potential should be less
% than this absolute value in normal operations).
rescaleY(h);
for ii=1:length(SCid)
  % Split YLabel to use two lines with "units" on the second line.
  ylabel(h(ii), strsplit(h(ii).YLabel.String,' '));
end;

if(bashRun)
  disp('Saving entire day');
  outPath = [dataPathRoot,'irfu',filesep,'plots',filesep];
  % Entire day:
 % savefig(fig,[outPath,DayOfInterest,'_000000_235959'],'compact'); %R2013b or later
  print(fig,'-dpng',[outPath,DayOfInterest,'_000000_235959']);
  % Interval 00:01:00 to 00:11:00 UTC
  tint = irf.tint([DayOfInterest,'T00:01:00Z'],[DayOfInterest,'T00:11:00Z']);
  irf_zoom(h,'x',tint);
  irf_zoom(h,'y'); rescaleY(h);
  title(h(1),{['Plot created: ',nowStr,'. DCE raw from our source cdf files for all four s/c.'],...
    'Zoom in: 00:01:00 to 00:11:00 UTC'});
  disp('Saving 00:01:00 to 00:11:00');
  print(fig,'-dpng',[outPath,DayOfInterest,'_000100_001100']);
  % Interval 06:00:00 to 06:10:00 UTC
  tint = irf.tint([DayOfInterest,'T06:00:00Z'],[DayOfInterest,'T06:10:00Z']);
  irf_zoom(h,'x',tint);
  irf_zoom(h,'y'); rescaleY(h);
  title(h(1),{['Plot created: ',nowStr,'. DCE raw from our source cdf files for all four s/c.'],...
    'Zoom in: 06:00:00 to 06:10:00 UTC.'});
  disp('Saving 06:00:00 to 06:10:00');
  print(fig,'-dpng',[outPath,DayOfInterest,'_060000_061000']);
  % Interval 12:00:00 to 12:01:00 UTC (ie a shorter interval as well).
  tint = irf.tint([DayOfInterest,'T12:00:00Z'],[DayOfInterest,'T12:01:00Z']);
  irf_zoom(h,'x',tint);
  irf_zoom(h,'y'); rescaleY(h);
  title(h(1),{['Plot created: ',nowStr,'. DCE raw from our source cdf files for all four s/c.'],...
    'Zoom in: 12:00:00 to 12:01:00 UTC.'});
  disp('Saving 12:00:00 to 12:01:00');
  print(fig,'-dpng',[outPath,DayOfInterest,'_120000_120100']);
  % Interval 18:00:00 to 18:10:00 UTC
  tint = irf.tint([DayOfInterest,'T18:00:00Z'],[DayOfInterest,'T18:10:00Z']);
  irf_zoom(h,'x',tint);
  irf_zoom(h,'y'); rescaleY(h);
  title(h(1),{['Plot created: ',nowStr,'. DCE raw from our source cdf files for all four s/c.'],...
    'Zoom in: 18:00:00 to 18:10:00 UTC.'});
  disp('Saving 18:00:00 to 18:10:00');
  print(fig,'-dpng',[outPath,DayOfInterest,'_180000_181000']);
end

  function rescaleY(h)
    narginchk(1,1);
    for id=1:length(SCid)
      % Get currently used limits
      ylimits = ylim(h(id));
      % If outside region -15 to +15, rescale.
      if(ylimits(1)<-15 || ylimits(2)>15)
        ylimits(1) = max(-15,ylimits(1));
        ylimits(2) = min(15, ylimits(2));
        ylim(h(id), ylimits);
      end;
    end
  end

end
