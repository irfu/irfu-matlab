% Engine for JUICE TM/power budget simulation.
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
% timeSec                : 1D array of time for which values should be evaluated.
% InsModes, RadModes     : 1D struct array with non-array fields (not counting strings):
%   .id                  : Formally defined string that identifies the mode.
%   .survDataBps         : Survey data production rate in bits/s, after compression.
%   .richDataBps         : Rich data production rate in bits/s, after compression.
%   .powerWatt           : Power consumption. Watt.
% InsModeSeq, RadModeSeq : 1D struct arrays with fields:
%   .id                  : ID string of instrument mode.
%   .beginSec            : Time the mode begins running.
% --
% NOTE: InsModeSeq/RadModeSeq(i).beginSec must increase in value with i.
% NOTE: All InsModes(i).id must be unique.
% NOTE: All RadModes(i).id must be unique.
% NOTE: Time format: Seconds from arbitrary epoch.
% 
%
% RETURN VALUE
% ============
% data : Non-array struct with 1D array fields.
%   .modeNbrIns      : Index into InsModes, representing the in situ mode.
%   .modeNbrRad      : Index into RadModes, representing the radio   mode.
%   .survDataBps
%   .survDataIntBits
%   .richDataBps
%   .richDataIntBits
%   .powerWatt
%
% 
% VARIABLE NAMING CONVENTIONS
% ===========================
% bps = Bits per second
% Int = Integrated over time
% Rad = Radio (mode)
% Ins = In situ (mode)
% Surv = Survey data (cf rich data)
% Seq = Sequence
%
%
% NOTES
% =====
% Not running any instrument mode is handled by defining and "running" a made-up "No-instrument-mode instrument mode".
% Only includes the power consumption from the in situ and radio modes, not from "powerDPU = 1912; powerXtra = 1100;" as
% in "ju_rpwi_tm_power.m".
%
%
% Created 2017-12-15 by Erik Johansson, IRF Uppsala.
%
function [data] = engine(timeSec, InsModes, RadModes, InsModeSeq, RadModeSeq)
% PROPOSAL: Change name of package.
%   PROPOSAL: Imply both TM and power.
%   PROPOSAL: Something with "budget". budget_simulator, budget_sim
%   PROPOSAL: Nothing with "simulator".
%   PROPOSAL: TM_power_budget, tm_power_budget
%   PROPOSAL: Something that implies JUICE?
%       PROPOSAL: Parent package named "juice"? JUICE? JUI?
% PROPOSAL: Change name of file/module/function.
%   PROPOSAL: Should imply "engine".
%
% TODO-DECISION: Use which terminology?
%   PROPOSAL: Types of data & TM: Survey data, rich data, stored survey data, stored rich data, 
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

% TODO-DECISION: Arguments
%       in situ modes, radio modes ("definitions")
%           production data rate, power consumption
%       List of when modes run.
%       On board storage: cluster size (not storage size)
%       Size of files produced.
%       Use as timestamps for return values?
% PROPOSAL: Return values:
%       Survey & rich data production rate, with and without accounting for file cluster size?
%       Power consumption as function of time.
%   
% PROPOSAL: Begin with ju_rpwi_tm_power.m, and rewrite it to use this engine.
% PROPOSAL: Automatic test code?


%===========================================================
% Create arrays of indices to in situ modes and radio modes
% ASSUMPTION: InsModeSeq(i).beginSec increases with i.
%===========================================================
[InsModeSeq.idNbr] = deal(0);   % Indices to mode.
[RadModeSeq.idNbr] = deal(0);   % Indices to mode.
for i = 1:numel(InsModeSeq)
    InsModeSeq(i).idNbr = find(strcmp(InsModeSeq(i).id, {InsModes.id}));
end
for i = 1:numel(RadModeSeq)
    RadModeSeq(i).idNbr = find(strcmp(RadModeSeq(i).id, {RadModes.id}));
end



[data.modeNbrIns, ~] = TM_power_budget.step_func([InsModeSeq.beginSec], [InsModeSeq.idNbr], timeSec);
[data.modeNbrRad, ~] = TM_power_budget.step_func([RadModeSeq.beginSec], [RadModeSeq.idNbr], timeSec);

[data.powerWatt, ~] = add_step_funcs(...
    [InsModeSeq.beginSec], [InsModes([InsModeSeq.idNbr]).powerWatt], ...
    [RadModeSeq.beginSec], [RadModes([RadModeSeq.idNbr]).powerWatt], timeSec);

[data.survDataBps, data.survDataIntBits] = add_step_funcs(...
    [InsModeSeq.beginSec], [InsModes([InsModeSeq.idNbr]).survDataBps], ...
    [RadModeSeq.beginSec], [RadModes([RadModeSeq.idNbr]).survDataBps], timeSec);

[data.richDataBps, data.richDataIntBits] = add_step_funcs(...
    [InsModeSeq.beginSec], [InsModes([InsModeSeq.idNbr]).richDataBps], ...
    [RadModeSeq.beginSec], [RadModes([RadModeSeq.idNbr]).richDataBps], timeSec);


%warning('UNFINISHED FUNCTION')
end



% Call TM_power_budget.step_func for different functions and add the results.
%
% Using this function simplifies and clarifies the main function a lot.
function [Y, YI] = add_step_funcs(xp1, yp1, xp2, yp2, X)
[Y1, YI1] = TM_power_budget.step_func(xp1, yp1, X);
[Y2, YI2] = TM_power_budget.step_func(xp2, yp2, X);
Y  = Y1  + Y2;
YI = YI1 + YI2;
end
