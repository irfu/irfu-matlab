function Plot_HK10E(DayOfInterest, bashRun)
% Plot_HK10E function for plotting bias settings (voltage and current) of
% all SDP probe elements on all spacecraft. Read from HK10E files.
%
% Note: will currently fail if HK10E files does not exist for all SC for
% the requested day.
%
%  Inputs:
%    DayOfInterest - day string in format 'YYYY-MM-DD'.
%    [bashRun]     - optional boolean to run in bash mode (save plots to
%           "$DATA_PATH_ROOT/irfu/plots" if it exist).
%  Usage:
%   Plot_HK10E('2015-05-12')
%
% Created: 2015/05/29
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
% SPDX-License-Identifier: CC-BY-NC-SA-4.0

narginchk(1,2);
if nargin == 1, bashRun=false; end
nowStr = char(datetime("now", "TimeZone","UTC", "Format","uuuu-MM-dd"));
if((numel(sscanf(DayOfInterest,'%4d-%2d-%2d%s'))~=3)||(length(DayOfInterest)~=10))
  error('Incorrect day string, format should be "YYYY-MM-DD".');
end

dataPathRoot = getenv('DATA_PATH_ROOT'); % Default to "/data/mms/"
dbList = mms_local_file_db(dataPathRoot);

if(nargin==2 && bashRun)
  outPath = ['/home/thoni/MMS',filesep,'plots',filesep];
  if(~exist(outPath,'dir')), error('outpath does not exist'); end
end

% HK_10E are daily files so start time is set to midnight of DayOfInterest
% and stop time midnight the next day.
tStart = [DayOfInterest,'T00:00:00.000000000Z'];
tStop  = [DayOfInterest,'T23:59:60.000000000Z'];
tint = irf.tint(tStart,tStop);

%% Identify and load all HK_10E files.
SCid = {'mms1', 'mms2', 'mms3', 'mms4'};
probes = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'};
for id=1:length(SCid)
  hk10eDB.(SCid{id}) = list_files(dbList,[SCid{id},'_fields_hk_l1b_10e'],tint);
  if(isempty(hk10eDB.(SCid{id})))
    warning(['No ',SCid{id},' HK_10E file found for this day.']);
    % Fill one datapoint with NaN to ensure plot function have something
    % to plot.
    for ii=1:length(probes)
      hk10eDB.(SCid{id}).dac.(probes{ii}) = irf.ts_scalar(tint.start, NaN);
      hk10eDB.(SCid{id}).og.(probes{ii}) = irf.ts_scalar(tint.start, NaN);
      hk10eDB.(SCid{id}).ig.(probes{ii}) = irf.ts_scalar(tint.start, NaN);
    end
  else
    % Load files
    hk10eDB.(SCid{id}).obj = dataobj([hk10eDB.(SCid{id}).path, filesep, hk10eDB.(SCid{id}).name]);
    toConvert = false;
    if tint.start < EpochTT('2015-06-01T00:00:00')
      % before this date some files are still old, (<0.5.z) with raw TM units
      fileVer = regexp(hk10eDB.(SCid{id}).name, ...
        'mms[1-4]_fields_hk_l1b_10e_\d{8,8}_v(?<verStr>\d{1,}.\d{1,}.\d{1,}).cdf', ...
        'names');
      if ~is_version_geq(fileVer.verStr, '0.5.0')
        % Old files, to be converted
        toConvert = true;
      end
    end
    for ii=1:length(probes)
      if(toConvert)
        % Convert them
        % DAC (probe current)
        hk10eDB.(SCid{id}).dac.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'dac']);
        hk10eDB.(SCid{id}).dac.(probes{ii}).data = dac2iBias(hk10eDB.(SCid{id}).dac.(probes{ii}).data);
        % OG (outer guard)
        hk10eDB.(SCid{id}).og.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'og']);
        hk10eDB.(SCid{id}).og.(probes{ii}).data = og2vBias(hk10eDB.(SCid{id}).og.(probes{ii}).data);
        % IG (inner guard)
        hk10eDB.(SCid{id}).ig.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'ig']);
        hk10eDB.(SCid{id}).ig.(probes{ii}).data = ig2vBias(hk10eDB.(SCid{id}).ig.(probes{ii}).data);
      else
        % Already converted.
        % DAC (probe current)
        hk10eDB.(SCid{id}).dac.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'dac']);
        % OG (outer guard)
        hk10eDB.(SCid{id}).og.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'og']);
        % IG (inner guard)
        hk10eDB.(SCid{id}).ig.(probes{ii}) = get_ts(hk10eDB.(SCid{id}).obj, [SCid{id}, '_10e_beb', num2str(ii), 'ig']);
      end
    end
  end
end

%% Plot All results.
% Note this will currently fail if not all HK10E files were located for
% all s/c for the day of interest.

% DAC
h = irf_plot(4, 'newfigure');
fig = gcf;
fig.WindowState = 'maximized';
for id=1:length(SCid)
  h(id) = irf_panel(SCid{id});
  irf_plot(h(id), {...
    hk10eDB.(SCid{id}).dac.P1, hk10eDB.(SCid{id}).dac.P2, ...
    hk10eDB.(SCid{id}).dac.P3, hk10eDB.(SCid{id}).dac.P4}, 'comp');
  legend(h(id), probes(1:4));
  ylabel(h(id), {[upper(SCid{id}), ' DAC'],'[nA]'});
  yMinMax = ylim(h(id));
  % adjust max/min by 6 to allow nicer looking plots (DAC has not ever
  % ended with xx6 so this will force plots to not align with panel lines.
  irf_zoom(h(id), 'y', [yMinMax(1)-6 yMinMax(2)+6]);
end
irf_zoom(h, 'x', tint);
title(h(1),['Plot created: ',nowStr,'. Probe current from HK\_10E for all four probes on all four s/c.']);

if(bashRun)
  disp('Saving HK 10E DAC.');
  pause(1);
  % print(fig, '-dpng',[outPath,DayOfInterest,'_10E_DAC']);
  exportgraphics(gcf, [outPath,DayOfInterest,'_10E_DAC.png']);
  % extract median DAB bias and save it
end

% OG
h = irf_plot(4, 'newfigure');
fig = gcf;
fig.WindowState = 'maximized';
for id=1:length(SCid)
  h(id) = irf_panel(SCid{id});
  irf_plot(h(id), {...
    hk10eDB.(SCid{id}).og.P1, hk10eDB.(SCid{id}).og.P2, ...
    hk10eDB.(SCid{id}).og.P3, hk10eDB.(SCid{id}).og.P4}, 'comp');
  legend(h(id), probes(1:4));
  ylabel(h(id), {[upper(SCid{id}), ' OG'],'[V]'});
  yMinMax = ylim(h(id));
  % adjust max/min by +/-0.1 to allow nicer looking plots.
  irf_zoom(h(id), 'y', [yMinMax(1)-0.1 yMinMax(2)+0.1]);
end
title(h(1),['Plot created: ',nowStr,'. Outer guard bias from HK\_10E for all four probes on all four s/c.']);
irf_zoom(h, 'x', tint);
if(bashRun)
  disp('Saving HK 10E OG.');
  pause(1);
  % print(gcf, '-dpng', [outPath,DayOfInterest,'_10E_OG']);
  exportgraphics(gcf, [outPath,DayOfInterest,'_10E_OG.png']);
end

% IG
h = irf_plot(4, 'newfigure');
fig = gcf;
fig.WindowState = 'maximized';
for id=1:length(SCid)
  h(id) = irf_panel(SCid{id});
  irf_plot(h(id), {...
    hk10eDB.(SCid{id}).ig.P1, hk10eDB.(SCid{id}).ig.P2, ...
    hk10eDB.(SCid{id}).ig.P3, hk10eDB.(SCid{id}).ig.P4}, 'comp');
  legend(h(id), probes(1:4));
  ylabel(h(id), {[upper(SCid{id}), ' IG'],'[V]'});
  yMinMax = ylim(h(id));
  % adjust max/min by +/-0.1 to allow nicer looking plots.
  irf_zoom(h(id), 'y', [yMinMax(1)-0.1 yMinMax(2)+0.1]);
end
title(h(1),['Plot created: ',nowStr,'. Inner guard bias from HK\_10E for all four probes on all four s/c.']);
irf_zoom(h, 'x', tint);
if(bashRun)
  disp('Saving HK 10E IG.');
  pause(1);
  exportgraphics(gcf, [outPath,DayOfInterest,'_10E_IG.png']);
  % print(gcf, '-dpng', [outPath,DayOfInterest,'_10E_IG']);
end

if(bashRun)
  % Save a median (per day) DAC value as a TSeries, to be compared with
  % analysed sweeps. Use median to avoid picking eclipse time when DAC is
  % zero or HK generated at the same time as sweep when DAC may be very
  % different from the nominal daily value.
  p1Dac=[]; p2Dac=[]; p3Dac=[]; p4Dac=[]; p5Dac=[]; p6Dac=[]; %#ok<NASGU>
  for id = 1:4
    c_eval('p?Dac = irf.ts_scalar(tint.start, median(hk10eDB.(SCid{id}).dac.P?.data));', 1:6);
    % do we have an existing "obj" file
    dacFile = fullfile('/data', 'mms', 'irfu', 'plots', 'edp', 'DAC', 'obj', [SCid{id}, '_dacTsCombined.mat']);
    if exist(dacFile, 'file')
      % load old values and combine TSeries
      load(dacFile, '-mat', ...
        'p1_dac', 'p2_dac','p3_dac','p4_dac', 'p5_dac', 'p6_dac');
      c_eval('p?_dac = combine(p?_dac, p?Dac);', 1:6);
    else
      % first run, just copy temp names to the combined name
      c_eval('p?_dac=p?Dac;', 1:6);
    end
    save(dacFile, ...
      'p1_dac', 'p2_dac', 'p3_dac', 'p4_dac', 'p5_dac', 'p6_dac', '-mat');
  end
end

%% Conversion functions, from TM to physical units.
% Source for conversion factors is "ADP/SDP BEB Command Interface", v1.7,
% 2012/02/02.
  function iBias = dac2iBias(dac)
    iBias = (double(dac)-32768)*25.28*10^-3; % nA
  end

  function vBias = og2vBias(og)
    vBias = (double(og)-32768)*317.5*10^-6; % V
  end

  function vBias = ig2vBias(ig)
    vBias = (double(ig)-32768)*317.5*10^-6; % V
  end
end
