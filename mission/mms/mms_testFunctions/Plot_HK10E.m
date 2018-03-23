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
%

narginchk(1,2);

nowStr = irf_time(now,'datenum>utc_yyyy-mm-dd');
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
probes = {'P1', 'P2', 'P3', 'P4'};
for id=1:length(SCid)
  hk10eDB.(SCid{id}) = list_files(dbList,[SCid{id},'_fields_hk_l1b_10e'],tint);
  if(isempty(hk10eDB.(SCid{id})))
    warning(['No ',SCid{id},' HK_10E file found for this day.']);
    % Fill one datapoint with NaN to ensure plot function have something
    % to plot.
    hk10eDB.(SCid{id}).time = tint.start.epochUnix;
    for ii=1:length(probes)
      hk10eDB.(SCid{id}).dac.(probes{ii}) = NaN;
      hk10eDB.(SCid{id}).og.(probes{ii}) = NaN;
      hk10eDB.(SCid{id}).ig.(probes{ii}) = NaN;
    end
  else
    % Load files
    hk10eDB.(SCid{id}).obj = dataobj([hk10eDB.(SCid{id}).path, filesep, hk10eDB.(SCid{id}).name]);
    % Convert to epoch times (for irf_plot)
    hk10eDB.(SCid{id}).time = irf_time(hk10eDB.(SCid{id}).obj.data.Epoch.data,'ttns>epoch');
    for ii=1:length(probes)
      if( tint.start < EpochTT('2015-06-01T00:00:00') )
        % Convert them
        % DAC (probe current)
        hk10eDB.(SCid{id}).dac.(probes{ii}) = dac2iBias(hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'dac']).data);
        % OG (outer guard)
        hk10eDB.(SCid{id}).og.(probes{ii}) = og2vBias(hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'og']).data);
        % IG (inner guard)
        hk10eDB.(SCid{id}).ig.(probes{ii}) = ig2vBias(hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'ig']).data);
      else
        % Already converted.
        % DAC (probe current)
        hk10eDB.(SCid{id}).dac.(probes{ii}) = hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'dac']).data;
        % OG (outer guard)
        hk10eDB.(SCid{id}).og.(probes{ii}) = hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'og']).data;
        % IG (inner guard)
        hk10eDB.(SCid{id}).ig.(probes{ii}) = hk10eDB.(SCid{id}).obj.data.([SCid{id},'_10e_beb',num2str(ii),'ig']).data;
      end
    end
  end
end

%% Plot All results.
% Note this will currently fail if not all HK10E files were located for
% all s/c for the day of interest.

% DAC
toPlot = cell(1,length(SCid));
for id=1:length(SCid)
  toPlot{1,id}=[hk10eDB.(SCid{id}).time, ...
    hk10eDB.(SCid{id}).dac.P1, ...
    hk10eDB.(SCid{id}).dac.P2, ...
    hk10eDB.(SCid{id}).dac.P3, ...
    hk10eDB.(SCid{id}).dac.P4];
end
fig=figure;
h = irf_plot(toPlot);
title(h(1),['Plot created: ',nowStr,'. Probe current from HK\_10E for all four probes on all four s/c.']);
for id=1:length(SCid)
  legend(h(id), probes);
  ylabel(h(id),{[upper(SCid{id}), ' DAC'],'[nA]'});
end
if(bashRun)
  disp('Saving HK 10E DAC.');
  print(fig,'-dpng',[outPath,DayOfInterest,'_10E_DAC']);
end

% OG
toPlot = cell(1,length(SCid));
for id=1:length(SCid)
  toPlot{1,id}=[hk10eDB.(SCid{id}).time, ...
    hk10eDB.(SCid{id}).og.P1, ...
    hk10eDB.(SCid{id}).og.P2, ...
    hk10eDB.(SCid{id}).og.P3, ...
    hk10eDB.(SCid{id}).og.P4];
end
fig=figure;
h = irf_plot(toPlot);
title(h(1),['Plot created: ',nowStr,'. Outer guard bias from HK\_10E for all four probes on all four s/c.']);
for id=1:length(SCid)
  legend(h(id), probes);
  ylabel(h(id),{[upper(SCid{id}), ' OG'],'[V]'});
end
if(bashRun)
  disp('Saving HK 10E OG.');
  print(fig,'-dpng',[outPath,DayOfInterest,'_10E_OG']);
end

% IG
toPlot = cell(1,length(SCid));
for id=1:length(SCid)
  toPlot{1,id}=[hk10eDB.(SCid{id}).time, ...
    hk10eDB.(SCid{id}).ig.P1, ...
    hk10eDB.(SCid{id}).ig.P2, ...
    hk10eDB.(SCid{id}).ig.P3, ...
    hk10eDB.(SCid{id}).ig.P4];
end
fig=figure;
h = irf_plot(toPlot);
title(h(1),['Plot created: ',nowStr,'. Inner guard bias from HK\_10E for all four probes on all four s/c.']);
for id=1:length(SCid)
  legend(h(id), probes);
  ylabel(h(id),{[upper(SCid{id}), ' IG'],'[V]'});
end
if(bashRun)
  disp('Saving HK 10E IG.');
  print(fig,'-dpng',[outPath,DayOfInterest,'_10E_IG']);
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
