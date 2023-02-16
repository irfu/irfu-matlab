% PLOT_DCE_V_PHASE reads and plots inital MMS data raw from source files.
%
%   This function was created in order to plot DCE, DCV and DEFATT Z-Phase
% from MMS in-flight data. It uses irfu-matlab and excerts from the MMS
% processing code, but not the MMS processing code itself as currently our
% data sources have issues with timestamps. (duplicated timestamps in
% DEFATT and going back and forth in time in DCE/DCV).
%
% Input (required):
%    DEFATT file
%    DCE cdf file
%    DCV cdf file
%
% Output (optional):
%    ZPhase,  Z-Phase from DefAtt file, [time (in epoch), ZPhase]
%    E1234,   Electric field from DCE file, [time (in epoch), E12, E34]
%    V1234,   Potential of each probe from DCV file, [time (in epoch), V1, V2, V3, V4]
%    dcePlot, handle for the DCE plot
%    dcvPlot, handle for the DCV plot
%
%   Example:
% [ZPHASE, E1234, V1234, dcePlot, dcvPlot] = Plot_DCE_V_Phase('/path/defattFile','/path/dceFile','/path/dcvFile');
%     or
% [ZPHASE, E1234, V1234] = Plot_DCE_V_Phase('/path/defattFile','/path/dceFile','/path/dcvFile');
%
% Created: 2015/03/26
% Author: T. Nilsson, IRFU
%
% License: CreativeCommons BY-NC-SA 4.0
% https://creativecommons.org/licenses/by-nc-sa/4.0/
% SPDX-License-Identifier: CC-BY-NC-SA-4.0

function [ZPHASE, E1234, V1234, dcePlot, dcvPlot] = Plot_DCE_V_Phase(defattFile, dceFile, dcvFile)

% Check input
narginchk(3,3);
global scID
scID = [];
checkFile(defattFile);
checkFile(dceFile);
checkFile(dcvFile);
% Check output
nargoutchk(0,5);

%% Read DEFATT
formatSpec = '%f-%f%s %*f %*f %*f %*f %*f %*f %*f %*f %f %*[^\n]';
headerGrep = 'COMMENT';
% Get number of last header line using unix commands grep, tail and cut.
if ismac
  [~, numHeaders] = unix(['grep -onr ',headerGrep,' ',defattFile,...
    ' | tail -n1 | cut -d'':'' -f2']);
else
  [~, numHeaders] = unix(['grep -onr ',headerGrep,' ',defattFile,...
    ' | tail -n1 | cut -d: -f1']);
end
numHeaders = str2double(numHeaders);

fileID = fopen(defattFile,'r');
tmpData = textscan( fileID, formatSpec,...
  'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', numHeaders );
fclose(fileID);
% Compute time
TTtime = [irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'), ...
  cell2mat(tmpData{1,3})];

% Store value to irf_plot [epoch, ZPhase]
ZPHASE = [irf_time(TTtime,'utc>epoch'), double(tmpData{1,4})];

clear TTtime fileID formatSpec tmpData numHeaders headerGrep

%% Read DCE file
dataObj = dataobj(dceFile);

% [epoch, E12, E34] for irf_plot
E1234 = [irf_time(dataObj.data.Epoch.data,'ttns>epoch'), ...
  double(dataObj.data.(['mms',scID,'_edp_dce_sensor']).data(:,1)), ...
  double(dataObj.data.(['mms',scID,'_edp_dce_sensor']).data(:,2))];

clear dataObj

%% Read DCV file
dataObj = dataobj(dcvFile);

% [epoch, V1, V2, V3, V4] for irf_plot
V1234 = [irf_time(dataObj.data.Epoch.data,'ttns>epoch'), ...
  double(dataObj.data.(['mms',scID,'_edp_dcv_sensor']).data(:,1)), ...
  double(dataObj.data.(['mms',scID,'_edp_dcv_sensor']).data(:,2)), ...
  double(dataObj.data.(['mms',scID,'_edp_dcv_sensor']).data(:,3)), ...
  double(dataObj.data.(['mms',scID,'_edp_dcv_sensor']).data(:,4))];

clear dataObj

%% Plot
today = char(datetime("now", "TimeZone","UTC", "Format","uuuu/MM/dd"));

% Set plot colors to same as that used for Cluster.
set(groot, 'defaultAxesColorOrder', [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]);

% DCE & ZPhase
figure; dcePlot = irf_plot({E1234, ZPHASE});
% Use strrep to replace "_" with "\_" in filename for title
title(dcePlot(1), {['Plot created ',today,', directly from cdf file'],...
  strrep(dceFile,'_','\_')});
ylabel(dcePlot(1), {'Electric field, Instrument ref. frame','[mV/m]'})
ylabel(dcePlot(2), {'Z-Phase, from DEFATT file','[deg]'})
legend(dcePlot(1), ['MMS',scID,' DCE12'],...
  ['MMS',scID,' DCE34']);
legend(dcePlot(2), 'Z-Phase (0 deg when sun is along BCS-X axis)');

% DCV & ZPhase
figure; dcvPlot = irf_plot({V1234, ZPHASE});
title(dcvPlot(1), {['Plot created ',today,', directly from cdf file'],...
  strrep(dcvFile,'_','\_')});
ylabel(dcvPlot(1), {'Probe potential, Ind. probes','[V]'})
ylabel(dcvPlot(2), {'Z-Phase, from DEFATT file','[deg]'})
legend(dcvPlot(1), ['MMS',scID,' DCV1'],...
  ['MMS',scID,' DCV2'],...
  ['MMS',scID,' DCV3'],...
  ['MMS',scID,' DCV4']);
legend(dcvPlot(2), 'Z-Phase (0 deg when sun is along BCS-X axis)');


%% Verify input file
  function checkFile(inFile)
    narginchk(1,1);
    if(exist(inFile,'file'))
      [~, scIDtmp, ~] = fileparts(inFile);
      scIDtmp = scIDtmp(4);
      if( ~isempty(scID) && ~strcmp(scID,scIDtmp) )
        error('Different SC');
      else
        scID = scIDtmp;
      end
    else
      error('File not found');
    end
  end

end
