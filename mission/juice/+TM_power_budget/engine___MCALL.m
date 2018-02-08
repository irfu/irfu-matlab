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

EngineConstants = TM_power_budget.default_constants();

tBegin = 0;
tEnd   = 24*3600;
switch 2
    case 1
        % GCO500. NOTE: Time in minutes. Later converted to seconds.
        seqInSitu = {
            0    'In-situ_slow';
            20   'In-situ_burst';
            25   'In-situ_slow';
            60   'In-situ_burst';
            65   'In-situ_slow';
            385  'In-situ_burst';
            390  'In-situ_slow';
            425  'In-situ_burst';
            430  'In-situ_slow';
            960  'In-situ_low';
            };
        seqRadio = {
            0    'Radio_full';
            62   'Radio_burst';
            72   'Radio_full';
            425  'Radar_mode-3';
            435  'Radio_full';
            500  'Radio_burst';
            510  'Radio_full';
            };
        seqClassif = {
            9*60, 2*2^20, 0;
            }; 
    case 2
        EngineConstants.InsModeDescrList(end+1) = struct('id', 'test_burst_1+1MiByps', 'prodSurvBps', 1*2^20*8, 'prodRichBps', 1*2^20*8, 'powerWatt', 5);
        seqInSitu = {
            0   'test_burst_1+1MiByps';
            };
        seqRadio = {
            0    'Radar_mode-3';
            };
        seqClassif = {
            6*60,   0*2^20,   0;
            };
        tEnd = 8*3600;  % Too high value triggers bug.
end
InsModeSeq = struct('beginSec', num2cell([seqInSitu{:, 1}]*60)', 'id', seqInSitu(:, 2));
RadModeSeq = struct('beginSec', num2cell([ seqRadio{:, 1}]*60)', 'id', seqRadio (:, 2));
ClassifSeq = struct('timeSec', num2cell([seqClassif{:, 1}]*60), 'selectedRichBytes', seqClassif(:, 2), 'rejectedRichBytes', seqClassif(:, 3));
tSec   = tBegin:10:tEnd;
tHours = tSec / 60^2;

CommEvents = [];
CommEvents.InsModeSeq  = InsModeSeq;
CommEvents.RadModeSeq  = RadModeSeq;
CommEvents.ClassifSeq  = ClassifSeq;
CommEvents.DownlinkSeq = struct('beginSec', {0}, 'bandwidthBps', {0});
InitialStorageState    = struct('queuedSurvBytes', 0, 'queuedRichBytes', 0, 'unclasRichBytes', 0);
[FinalStorageState, StateArrays, Clf] = TM_power_budget.engine(EngineConstants, InitialStorageState, CommEvents, tBegin, tEnd, tSec);

close all
XTICK = [0:3:24];
data = StateArrays;



if 1
    %=====
    figure
    %=====
    h = subplot(2,1,1);
    plot(h, tHours, data.iInsModeDescr, '.');
    ylabel(h, 'In situ mode')
    set(h, 'XTick', XTICK)
    set(h, 'XTickLabel', [])
    set(h, 'YTick', [min(data.iInsModeDescr):max(data.iInsModeDescr)])
    
    h = subplot(2,1,2);
    plot(h, tHours, data.iRadModeDescr, '.');
    xlabel('Time [h]')
    ylabel(h, 'Radio mode')
    set(h, 'XTick', XTICK)
    set(h, 'YTick', [min(data.iRadModeDescr):max(data.iRadModeDescr)])
    
    print('-dpng', 'Fig_modes.png')
end

if 1
    %=====
    figure
    %=====
    h = subplot(3,1,1);
    %semilogy(h, tHours, data.prodSurvBps+data.prodRichBps)
    semilogy(h, tHours, [data.prodSurvBps', data.prodRichBps'], '-'); legend('Survey', 'Rich')
    ylabel(h, 'Data production [bits/s]')
    set(h, 'XTickLabel', [])
    set(h, 'XTick', XTICK)
    
    h = subplot(3,1,2);
    usedStorageMiB = data.usedStorageBytes / 2^20;
    storageUnclasRichMiB = data.unclasRichBytes / 2^20 * Clf.rich;
    storageQueuedRichMiB = data.queuedRichBytes / 2^20 * Clf.rich;
    storageQueuedSurvMiB = data.queuedSurvBytes / 2^20 * Clf.surv;
    plot(h, tHours, [usedStorageMiB; storageUnclasRichMiB; storageQueuedRichMiB; storageQueuedSurvMiB])
    legend('Total', 'Unclas. rich', 'Queued rich', 'queued survey')
    %plot(h, tHours, [data.survDataIntBits', data.richDataIntBits']/2^20/8); legend('Survey', 'Rich')
    ylabel(h, 'Used storage [MiB]')
    set(h, 'XTickLabel', [])
    set(h, 'XTick', XTICK)
    
    h = subplot(3,1,3);
    plot(h, tHours, data.powerWatt)
    xlabel('Time [h]')
    ylabel('Power [W]')
    set(h, 'XTick', XTICK)
    
    print('-dpng', 'Fig_TM.png')
end

end
