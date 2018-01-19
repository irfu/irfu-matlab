% "Manual" test/demo of using other files in this package. MCALL = Manual call (convention).
%
%
%
% Created 2017-12-15 by Erik Johansson, IRF Uppsala.
%
function data = engine___MCALL
% PROPOSAL: Replace file by proper "demo"?
% PROPOSAL: Replace file by proper "scenario"?
% PROPOSAL: Add title to figures, ~"GCO500".
%=======================================================================================
% Copied from Dropbox/JUICE/Engineering/Budget_Tables/ju_rpwi_tm_power.m 2018-01-02.
%   case 'GCO500'
%     seqIns = [
%       0    2; ...
%       20   4; ...
%       25   2; ...
%       60   4; ...
%       65   2; ...
%       385  4; ...
%       390  2; ...
%       425  4;
%       430  2;
%       960  1; ...
%       ];
%
%     seqRad = [
%       0    6; ...
%       62   7; ...
%       72   6; ...
%       425  8; ...
%       435 6; ...
%       500  7; ...
%       510  6; ...
%       ];
%     case 'GCO500ext_GC0200'
%       modeTitle = 'GCO500ext/GC0200';
%
%       seqIns = [
%         0    3; ...
%         20   4; ...
%         25   3; ...
%         425  4;
%         430  3;
%         960  1; ...
%         ];
%
%       seqRad = [
%         0    6; ...
%         425  8; ...
%         435 6; ...
%         ];
%=======================================================================================

[InSituModes, RadioModes] = TM_power_budget.default_instrument_modes();

% GCO500. NOTE: Time in minutes. Later converted to seconds.
seqInSitu = {
    0    'In-situ_slow'; ...   % 5.489 W
    20   'In-situ_burst'; ...  % 5.489 W
    25   'In-situ_slow'; ...
    60   'In-situ_burst'; ...
    65   'In-situ_slow'; ...
    385  'In-situ_burst'; ...
    390  'In-situ_slow'; ...
    425  'In-situ_burst';
    430  'In-situ_slow';
    960  'In-situ_low'; ...    % 
    };
seqRadio = {
    0    'Radio_full'; ...   % 1.673 W
    62   'Radio_burst'; ...  % 6.149 W
    72   'Radio_full'; ...
    425  'Radar_mode-3'; ... % 
    435  'Radio_full'; ...
    500  'Radio_burst'; ...
    510  'Radio_full'; ...
    };
InSituModeSequence = struct('beginSec', num2cell([seqInSitu{:, 1}]*60)', 'id', seqInSitu(:, 2));
RadioModeSequence  = struct('beginSec', num2cell([ seqRadio{:, 1}]*60)', 'id', seqRadio (:, 2));

tSec   = 0:60:24*60*60;
tHours = tSec / 60^2;

data = TM_power_budget.engine(tSec, InSituModes, RadioModes, InSituModeSequence, RadioModeSequence);

close all
XTICK = [0 3 6 9 12 15 18 21 24];



%=====
figure
%=====
h = subplot(2,1,1);
plot(h, tHours, data.modeNbrIns, '.');
ylabel(h, 'In situ mode')
set(h, 'XTick', XTICK)
set(h, 'XTickLabel', [])
set(h, 'YTick', [min(data.modeNbrIns):max(data.modeNbrIns)])

h = subplot(2,1,2);
plot(h, tHours, data.modeNbrRad, '.');
ylabel(h, 'Radio mode')
set(h, 'XTick', XTICK)
set(h, 'YTick', [min(data.modeNbrRad):max(data.modeNbrRad)])

print('-dpng', 'Fig_modes.png')



%=====
figure
%=====
h = subplot(3,1,1);
%semilogy(h, tHours, data.survDataBps+data.richDataBps)
semilogy(h, tHours, [data.survDataBps', data.richDataBps'], '-'); legend('Survey', 'Rich')
ylabel(h, 'TM [bits/s]')
set(h, 'XTickLabel', [])
set(h, 'XTick', XTICK)

h = subplot(3,1,2);
%plot(h, tHours, (data.survDataIntBits+data.richDataIntBits)/2^20/8)
plot(h, tHours, [data.survDataIntBits', data.richDataIntBits']/2^20/8); legend('Survey', 'Rich')
ylabel(h, 'Integrated TM [MiB]')
set(h, 'XTickLabel', [])
set(h, 'XTick', XTICK)

h = subplot(3,1,3);
plot(h, tHours, data.powerWatt)
xlabel('Time [h]')
ylabel('Power [W]')
set(h, 'XTick', XTICK)

print('-dpng', 'Fig_TM.png')


end
