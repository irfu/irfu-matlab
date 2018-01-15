%
% Define constants representing properties for instrument modes needed when estimating TM and power use.
% The format is defined by TM_power_budget.engine.
%
% This function is meant to collect definitions of frequently used hardcoded values and do nothing else.
%
%
% Created 2018-01-02 by Erik Johansson, IRF Uppsala.
%
function [InSituModes, RadioModes] = default_instrument_modes
% TODO: Separate survey data and rich data.
% Survey and rich data production rates from
% https://docs.google.com/spreadsheets/d/1OTSdC_eI7N-mA29RaMZR-VnD4CaauRnuR8f4rVGBBoY/edit#gid=1675236151
% on the assumption that "in_situ-slow" contains only and all (in situ mode) survey data.
%
% Copied from ju_rpwi_tm_power.m 2018-01-02. Power values have been taken from there.
% ----------------------------------------------------------------
% Modes = {'In-situ_low',... %1
%   'In-situ_slow',...
%   'In-situ_normal',...
%   'In-situ_burst',...
%   'Radio_low',... % 5
%   'Radio_full',...
%   'Radio_burst',...
%   'Radar_mode-3'};
% bpsSC =   [ 979, 2097, 9806,65483,  135,  230,43200,143616];
% powerSC = [2559, 5489, 5489, 5489, 1673, 1673, 6149,  6149];
% ----------------------------------------------------------------

EmptyStructArray = struct('id', {}, 'survDataBps', {}, 'richDataBps', {}, 'powerWatt', {});

InSituModes = EmptyStructArray;
InSituModes(end+1) = create_record('In-situ_low',     979, 0, 2.559);
InSituModes(end+1) = create_record('In-situ_slow',   2097, 0, 5.489);
InSituModes(end+1) = create_record('In-situ_normal', (90+112+120+170), (1669+ 3823+3058+ 765), 5.489);
InSituModes(end+1) = create_record('In-situ_burst',  (90+112+480+170), (1669+53406+7646+1911), 5.489);

RadioModes = EmptyStructArray;
RadioModes(end+1) = create_record('Radio_low',       135, 0, 1.673);
RadioModes(end+1) = create_record('Radio_full',      230, 0, 1.673);
RadioModes(end+1) = create_record('Radio_burst',   43200, 0, 6.149);
RadioModes(end+1) = create_record('Radar_mode-3', 143616, 0, 6.149);
end



function m = create_record(id, survDataBps, richDataBps, powerWatt)
m.id          = id;
m.survDataBps = survDataBps;
m.richDataBps = richDataBps;
m.powerWatt   = powerWatt;
end
