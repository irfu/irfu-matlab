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
% Copied from ju_rpwi_tm_power.m 2018-01-02. Values have been taken from there.
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

EmptyStructArray = struct('id', {}, 'prodRateBps', {}, 'powerConsWatt', {});

InSituModes = EmptyStructArray;
InSituModes(end+1) = create_record('In-situ_low',     979, 2559);
InSituModes(end+1) = create_record('In-situ_slow',   2097, 5489);
InSituModes(end+1) = create_record('In-situ_normal', 9806, 5489);
InSituModes(end+1) = create_record('In-situ_burst', 65483, 5489);

RadioModes = EmptyStructArray;
RadioModes(end+1)   = create_record('Radio_low',       135, 1673);
RadioModes(end+1)   = create_record('Radio_full',      230, 1673);
RadioModes(end+1)   = create_record('Radio_burst',   43200, 6149);
RadioModes(end+1)   = create_record('Radar_mode-3', 143616, 6149);
end



function m = create_record(id, prodRateBps, powerConsWatt)
m.id = id;
m.prodRateBps = prodRateBps;
m.powerConsWatt = powerConsWatt;
end
