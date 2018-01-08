% Engine in JUICE TM/POWER/BUDGET simulator.
%
% UNFINISHED
%
%
% INTENTION FOR CODE
% ==================
% The main purpose for this code is to function as an "engine" for estimating/calculating:
% (1) TM (survey data produced, rich data produced, storage used, amount of rich data that can be downlinked), and
% (2) (possibly) power use over time.
% The caller supplies definitions of instrument modes and time series of when modes run.
%
% The code is meant to be an "engine", i.e. only solve an abstract task and not have any human-user-friendly UI, in
% particular no GUI. Such UIs should be provided by other code that calls this code.
%
% The engine is not meant to contain any hardcoded constants or default values.
% The engine is not meant to warn of exceeding budgets etc since (1) that is easy for calling code to do itself, and (2)
% may need to be done in a custom way.
%
% The code is meant to begin simple, and become more advanced over time.
% Therefore, no guarantuees of backward compatibility can be made for a long time.
% A future version can hopefully be used by official future applications.
%
% 
% ARGUMENTS
% =========
% InSituModes, RadioModes : 1D struct arrays with fields:
%   .id            : Formally defined string that identifies the mode.
%   .prodRateBps   : Production data rate in bits/s, after compression.
%   .powerConsWatt : Power consumption. Watt.
% InSituModeSequence, RadioModeSequence : 1D struct arrays with fields:
%   .beginSec : Time the mode begins running. Seconds from arbitrary epoch.
%   .id       : ID string of instrument mode.
% 
%
% TERMINOLOGY
% ===========
% bps = bits per second
%
%
% Created 2017-12-15 by Erik Johansson, IRF Uppsala.
%
function [] = engine(InSituModes, RadioModes, InSituModeSequence, RadioModeSequence)
% PROPOSAL: Change name of package.
%   PROPOSAL: Imply both TM and power.
%   PROPOSAL: Something with "budget". budget_simulator, budget_sim
%   PROPOSAL: Nothing with "simulator".
%   PROPOSAL: TM_power_budget, tm_power_budget
% PROPOSAL: Change name of file/module/function.
%   PROPOSAL: Should imply "engine".
%
% TODO-DECISION: Use which time format?
%
% TODO-DECISION: Use which terminology?
%   PROPOSAL: Types of data & TM: Survey data, rich data, stored survey data, stored rich data, 
%   PROPOSAL: Types of modes: ins=in situ, rad=radio
%       TODO-DECISION: Camel case for "inSitu" in identifiers?
%   PROPOSAL: Define term for ~"instrument mode". Mode?
%   PROPOSAL: Instrument mode ID strings
%   PROPOSAL: Standard abbreviations to be used in identifiers.
%
% TODO-DECISION: Define (and document) exact model for TM and power.
%   TODO-NEED-INFO: How/when is stored rich data deleted?
%   TODO-NEED-INFO: Only one power budget? (one scalar?)
%   TODO-NEED-INFO: Is no input variable truly ever a response to an output variable?
%       Can the problem be modelled as function: input-->output?
%       Ex: Amount of rich data downlinked.
%       Ex: Can not "derive" which modes to use depending on available downlink, power.
%   TODO-DECISION/NEED-INFO: Accuracy of model?
%       How handle time in return values?
%           (1) Need exact times in output (begin-end of "behaviours"/modes), or
%           (2) enough to have specific time resolution (state at specific times)?
%   TODO-DECISION: How handle compression?
%       PROPOSAL: Only work with uncompressed data.
%       PROPOSAL: Specify compression ratio. Different for each mode?
%           NOTE: Different compression ratios for different data products can be converted into one compression ratio
%                 for all data produced by a specific instrument mode.
%
% TODO-DECISION: Use 1D struct arrays, or structs with array fields?
%
% NEED: Ability to not run any instrument mode.
%   PROPOSAL: Define instrument mode that represents this.
%
% PROPOSAL: Arguments
%       in situ modes, radio modes ("definitions")
%           production data rate, power consumption
%       List of when modes run.
%       On board storage: cluster size (not storage size)
%       Size of files produced.
%   PROPOSAL: No data for which to compare results with
%       Ex: Available power over time
%       Ex: Available downlink TM over time.
%       Ex: Available onboard storage size.
% PROPOSAL: Return values: Survey/rich data production rate, power consumption as functions of time.
%   PROPOSAL: Ackumulated (integrated) production rates over time.
%   
% PROPOSAL: Begin with ju_rpwi_tm_power.m, and rewrite it to use this engine.
% PROPOSAL: (at least initially) define separe function with hardcoded values for instrument modes on format which can
%           easily be used for this function.


warning('UNFINISHED FUNCTION')
end
