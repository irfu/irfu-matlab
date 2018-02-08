function [StorageState1, StateArrays, Clf] = engine(Constants, StorageState0, CommEvents, timeSec0, timeSec1, timeArraySec)
% Engine for JUICE TM/power budget simulation.
%
%
% INTENTION FOR CODE
% ==================
% The main purpose of this code is to function as an "engine" for estimating/calculating/simulating:
% (1) The "flow of data" i.e.
%   (a) amount of data produced when
%   (b) how much data is stored onboard, 
%   (c) how long data is be stored onboard before it has to be deleted (rejected or downlinked)
%   (d) how much data is be downlinked
% and, possibly
% (2) power use over time
% given user-supplied
% (1) initial state
% (2) sequences of commanded "events" over time: beginning instrument modes, selection & rejection of rich data,
% downlink bandwidth.
%
%
% IMPLEMENTATION NOTES
% ====================
% * The code is meant to be an "engine" only, i.e. only solve an abstract task and not have any human-user-friendly UI,
% in particular NOT a GUI but also NOT a UI for easy use from the MATLAB command-line. Such UIs should be provided by
% other code that calls this code.
% * The engine is not meant to contain any hardcoded constants or default values.
% * The code is designed so that the storage state can be evolved in steps, via multiple chained calls to this function,
% i.e. from t_1 to t_2, then from t_2 to t_3 and so on.
% * The code is meant to begin as a simple version, and become more advanced over time. Therefore, no guarantuees of
% backward compatibility can be made for a long time. A future version can hopefully be used by official future
% applications.
%
%
% ARGUMENTS
% =========
% Constants                                 : Struct
%   .InsModeDescrList, .RadModeDescrList    : 1D struct array with non-array fields (not counting strings):
%       .id                                 : Formally defined string that identifies the mode.
%       .prodSurvBps                        : Survey data production rate. bits/s, after compression.
%       .prodRichBps                        : Rich   data production rate. bits/s, after compression.
%       .powerWatt                          : Power consumption.
%   .SystemPrps                             : Struct with constants definining the static properties of the system.
%       .storageBytes                : Total size of onboard storage memory available for storing instrument data (rich & survey). 
%       .clusterSizeBytes                   : Size of which every file is an integer multiple.
%       .survDataFileSizeBytes              : Influences cluster loss.
%       .richDataFileSizeBytes              : Influences cluster loss.
%
% InitialStorageState                       : State of the storage at time timeSec0.
%   .queuedSurvBytes                        : Stored survey data that is waiting to be downlinked (automatically all survey data).
%   .queuedRichBytes                        : Stored rich   data that is waiting to be downlinked (and has already been
%                                             selected for being downlinked).
%   .unclasRichBytes                        : Stored rich   data that is waiting for being selected/rejected for downlink.
%                                             (unclas=unclassified).
%                                             NOTE: All these measures refer to amount of data stored, NOT the amount of
%                                             storage space used (which is subject to cluster loss).
%
% CommEvents                                : Struct (not 1D array) describing (tele)commanded events occurring at specified timestamps.
%   .InsModeSeq, .RadModeSeq                : 1D struct arrays.
%       .beginSec                           : Time the mode begins running. Must increase monotonically.
%       .id                                 : ID string of instrument mode. Must be unique (separately for
%                                             InsModeSeq/RadModeSeq).
%   .ClassifSeq                             : 1D struct array. Each element describes an instant event consisting of unclassified
%                                             rich data being rejected (deleted) and/or selected (queued for downlink).
%       .timeSec                            : Time stamp for event.
%       .selectedRichBytes
%       .rejectedRichBytes
%   .DownlinkSeq                            : 1D struct array
%       .beginSec                           : Time after which the downlink is set to a certain constant bandwidth.
%       .bandwidthBps                       : Downlink bandwidth. Bits/s. Assuming no additional compression for the downlink.
%
% timeSec0, timeSec1                        : Start and end time for evolving the state.
% timeArraySec                              : 1D array of time for which values should be evaluated.
% --
% NOTE: All InsModeDescrList(i).id must be unique.
% NOTE: All RadModeDescrList(i).id must be unique.
% NOTE: Time format for all time stamps: Seconds ("Sec") from one common arbitrary epoch. (All time stamps are
% spacecraft time.)
% 
%
% RETURN VALUES
% =============
% StorageState1         : Same format as StorageState0. Valid at timeSec0.
% StateArrays           : Non-array struct with 1D array fields. Each array element is valid at time timeArraySec.
%   .timeArraySec
%   .iInsModeDescr      : Index into InsModeDescrList, representing the in situ mode.
%   .iRadModeDescr      : Index into RadModeDescrList, representing the radio   mode.
%   .unclasRichBytes
%   .queuedRichBytes
%   .queuedSurvBytes
%   .prodSurvBps
%   .prodRichBps
%   .usedStorageBytes   : Not on StorageStates, since can be derived from it (redundant).
%   .powerWatt
% Clf                   : Struct with cluster loss factors (>=1), the ratio between allocated storage and data in file.
%   .surv
%   .rich
%
%
% VARIABLE NAMING CONVENTIONS
% ===========================
% bps = Bits per second (not BYTES per second)
% Int = Integrated over time
% Rad = Radio (mode)
% Ins = In situ (mode)
% Surv = Survey (data)
% Rich = Rich (data)
% Seq = Sequence
%
%
% NOTES
% =====
% * Not running any instrument mode is modelled by defining and "running" a made-up "No-instrument-mode instrument mode"
% which produces zero data.
% * Not having any downlink is modelled by setting Downlink(i).bandwidthBps == 0.
% * Only includes the power consumption from the in situ and radio modes, not from "powerDPU = 1912; powerXtra = 1100;" as
% in "ju_rpwi_tm_power.m".
%
%
% IMPLEMENTATION NOTES
% ====================
% * Amount of stored data is specified in bytes since MB, GB etc would be slightly ambiguous.
% * Only permitting two file sizes, for survey and rich data separately, to keep the model simple. Having more file sizes
% (e.g. one per mode) complicates the representation of stored data since the relationship between (1) amount of
% data stored, and (2) the amount of storage space used, is no longer a constant ratio, i.e. one has to separately
% represent data with different file sizes with different scalars. Additionally, it then becomes important what is the
% file size of the data that is specified in the TC when it implies deleting data from onboard storage (downlink, and
% when selecting rich data), so that one can compute the amount of freed onboard storage space.
% 
%
% Initially created 2017-12-15 by Erik Johansson, IRF Uppsala.
%

% PROPOSAL: Change name of package.
%   PROPOSAL: Imply both TM and power.
%   PROPOSAL: Something with "budget". budget_simulator, budget_sim
%   PROPOSAL: Something that implies evolving over time: simulator, evolve, update.
%   PROPOSAL: TM_power_budget, tm_power_budget
%   PROPOSAL: Something that implies JUICE.
%       PROPOSAL: Parent package named "juice"? JUICE? JUI?
% PROPOSAL: Change name of file/module/function.
%   PROPOSAL: Should imply "engine".
%
% TODO-DECISION: Define (and document) exact model for TM and power.
%   TODO-NEED-INFO: Only one power budget? (one scalar?)
%   TODO-NEED-INFO: Is no input variable truly ever a response to an output variable?
%       Can the problem be modelled as function: input-->output?
%       Ex: Amount of rich data downlinked.
%       Ex: Can not "derive" which modes to use depending on available downlink, power.
%   TODO-DECISION: How handle compression?
%       PROPOSAL: Only work with uncompressed data.
%       PROPOSAL: Specify compression ratio. Different for each mode?
%           NOTE: Different compression ratios for different data products can be converted into one compression ratio
%                 for all data produced by a specific instrument mode.
%       PROPOSAL: Have upper & lower limit depending on compression.
%
% PROPOSAL: Add for every mode: File size for survey & rich data respectively (two numbers).
%
% TODO-DECISION: Use 1D struct arrays, or structs with array fields?
%   NOTE: Can not (conveniently) use 1D struct arrays for structs within structs. (?)
%
% PROPOSAL: Begin with ju_rpwi_tm_power.m, and rewrite it to use this engine.
% PROPOSAL: Automatic test code?
%
% PROPOSAL: Somehow add warning for "error states" (term?).
%
% NOTE: There are discrete events that can occur which are not in the preprogrammed sequences.
%   NOTE: Non-error events. ==> Must be able to handle to avoid negative queued data values.
%       (1) Runs out of queued survey data to download. ==> Use excess downlink to downlink selected rich data
%   NOTE: Error states
%       (1) Runs out of queued survey & rich data to download. ==> Excess downlink unused.
%       (2) Runs out of storage space. ==> Stop filling up storage space (rich data first, then survey?)
%
% PROPOSAL: Change to an implementation model which works in steps of time sequences uninterrupted by "events", i.e. time
%           sequences, i.e. mode changes, downlink changes, select/reject rich data events.
%   PRO: Can more naturally introduce new discrete "events".
%   ~CON: Must merge together the information from multiple sequences to produce long arrays of ~state values.
%   TODO-DECISION: How create long, multi-interval state value arrays.
%       PROPOSAL: Interpolate.
%           NOTE: Different kinds of interpolation for different data: Linear interpolation, latest-value interpolation.
%           PRO: Can be done once, after end of interval loop using all the storage states (at interval ends).
%
% PROPOSAL: Use numerical methods to solve differential equations.
%   TODO-DECISION: How generalize to multiple file sizes? Arbitrary files size changes? Arbitrary rules for downlinking
%   different file size data?
%   CON: Numerical methods do not work for discontinuous functions?
%
% PROPOSAL: Use specially crafted functions for working with step functions (and maybe ramp functions, min-max function
%           etc) to solve differential equations analytically.
%   CON: Can not solve for M_s (stored survey data) etc since occurs twice in equation, of which one occurrence is derivative.
%
% PROPOSAL: Have the caller set cluster loss factors, INSTEAD OF files sizes and cluster size.
%
% PROPOSAL: Return quantities for creating "XB-like corridor plot".
%   (1) Lower limit = Accumulated downlink integrated.
%   (2) Upper limit = Lower limit + storage size
%   (3) Data = Ackumulated data production
%   TODO-DECISION: Unclear how to handle rich data classification.
%       PROPOSAL: Let ackumulated data decrease at data rejection.
%       PROPOSAL: Let ackumulated downlink include rejected data.
%           CON: Not as predictable as downlink sequence.
%
% PROPOSAL: Ignore error states for now.
%   CON/NOTE: Still not enough to only rely on step_function because of
%       (1) data rejection. ==> ~downlink-like spikes
%       (2) Queud survey data > 0 still creates discrete events (beyond time sequence events)
%
% PROPOSAL: Only use bytes per second.
%   PRO: Consistent
%   PRO: Less risk of bugs.
% PROPOSAL: Shorten notation: Q=Queued, U=Unclassified, ClusterLossFactor-->aboveMax
% TODO: Clarify how handling commanded events at beginning and end. Needed for chaining.


% Variable naming conventions
% ---------------------------
% S = State
% (0,1 = Initial and final state/time for the call to this function)
% 2,3 = Refers to beginning and end of time interval which includes at most one commanded event which can only take place at time 3.
% CLF = Cluster loss factor (>=1)

%========================
% Shorten variable names
%========================
SystemPrps       = Constants.SystemPrps;
InsModeDescrList = Constants.InsModeDescrList;
RadModeDescrList = Constants.RadModeDescrList;
clear Constants;
InsModeSeq  = CommEvents.InsModeSeq;
RadModeSeq  = CommEvents.RadModeSeq;
DownlinkSeq = CommEvents.DownlinkSeq;
ClassifSeq  = CommEvents.ClassifSeq;
clear CommEvents;


EMPTY_ARRAY_FIELDS = {'unclasRichBytes', 'queuedRichBytes', 'queuedSurvBytes', 'iInsModeDescr', 'iRadModeDescr', 'prodSurvBps', 'prodRichBps', 'powerWatt', 'usedStorageBytes', 'downlinkBps'};
tempEmptyArray = zeros(size(timeArraySec)) * NaN;
StateArrays = [];
StateArrays.timeArraySec = timeArraySec;
for i = 1:numel(EMPTY_ARRAY_FIELDS)
    StateArrays.(EMPTY_ARRAY_FIELDS{i}) = tempEmptyArray;
end

SystemPrps.survClf = cluster_loss_factor(SystemPrps.survFileSizeBytes, SystemPrps.clusterSizeBytes);
SystemPrps.richClf = cluster_loss_factor(SystemPrps.richFileSizeBytes, SystemPrps.clusterSizeBytes);

% Create one single, time-ordered list of commanded events.
[commEventTimesSec, iEventType, jEvent] = merge_seq([InsModeSeq.beginSec], [RadModeSeq.beginSec], [ClassifSeq.timeSec], [DownlinkSeq.beginSec]);


% Find first commanded event.
iCommEvent = find(commEventTimesSec >= timeSec0, true, 'first');



%=======================================================
% Evolve state 0-->1
% Split evolution into smaller steps 2-->3.
%=======================================================
S2 = StorageState0;
timeSec2 = timeSec0;

% Find initial mode and downlink
% NOTE: It does not matter if setting in situ/radio mode or downlink twice (at the same timestamp) as opposed to
% classification events.
S2.iInsModeSeq  = prev_lower_equal_value([InsModeSeq.beginSec],  timeSec0);
S2.iRadModeSeq  = prev_lower_equal_value([RadModeSeq.beginSec],  timeSec0);
S2.iDownlinkSeq = prev_lower_equal_value([DownlinkSeq.beginSec], timeSec0);
% ASSERTIONS
if isempty(S2.iInsModeSeq)
    error('Can not find in-situ mode at beginning of time interval.')
elseif isempty(S2.iRadModeSeq)
    error('Can not find radio mode at beginning of time interval.')
elseif isempty(S2.iDownlinkSeq)
    error('Can not find downlink at beginning of time interval.')
end

% Obtain actual mode descriptions and downlink value.
% TODO-NEED-INFO: Are these modes descriptions at all really used??!!!
[S2.InsModeDescr, S2.iInsModeDescr] = get_mode_descr(InsModeDescrList, InsModeSeq(S2.iInsModeSeq).id);
[S2.RadModeDescr, S2.iRadModeDescr] = get_mode_descr(RadModeDescrList, RadModeSeq(S2.iRadModeSeq).id);
S2.downlinkBps = DownlinkSeq(S2.iDownlinkSeq).bandwidthBps;

lastTimeInterval = false;
while true
    % BUG?!!: Need way of exiting while loop if timeSec3 == timeSec1?!!

    % Set end of interval (timestamp 3), generally being the next commanded event.
    if (iCommEvent > numel(commEventTimesSec)) || (commEventTimesSec(iCommEvent) > timeSec1)
        % CASE: (1) There is no more commanded event, or (2) the next commanded event takes place after timeSec1.
        timeSec3         = timeSec1;
        lastTimeInterval = true;
    else
        timeSec3 = commEventTimesSec(iCommEvent);
    end

    %===========================================
    % Evolve system 2-->3
    % Split evolution into smaller steps 4-->5.
    %===========================================
    timeSec4 = timeSec2;
    timeSec5 = timeSec3;
    S4 = S2;
    while true
        %S4 = complement_state(S4, SystemPrps);
    
        %=============================================
        % Try to evolve system into the future, 4-->5
        %=============================================
        %timeSec4
        %timeSec5
        [S5, timeSecP, StateArrays] = evolve_state_linearly(SystemPrps, S4, timeSec4, timeSec5, StateArrays);
        
        % DEBUG
        %S5 = complement_state(S5, SystemPrps);
        % Find time C, where system reached some sort of "bounds".
        %timeSecC = timeSec5;
            
        % Decide what to do next.
        if timeSecP < timeSec5
            % CASE: Did NOT evolve all the way.
            timeSec5 = timeSecP;
            % ==> Try evolving again, but for shorter time interval.
        else
            % CASE: Did succeed in evolving until timeSec5.
            if timeSec3 == timeSec5
                % CASE: Did evolve until timeSec3.
                S3 = S5;
                break
            else
                % CASE: Did not evolve until timeSec3.
                S4 = S5;
                timeSec4 = timeSec5;
                timeSec5 = timeSec3;
                % ==> Try evolving again, but starting at 4.
            end
        end
    end    
    
    if lastTimeInterval
        break
    end
    
    %===========================================================
    % Handle event specified in caller-submitted time sequences
    %===========================================================
    j3 = jEvent(iCommEvent);
    switch iEventType(iCommEvent)
        
        case 1  % In situ mode change
            
            id = InsModeSeq(j3).id;
            [S3.InsModeDescr, S3.iInsModeDescr] = get_mode_descr(InsModeDescrList, id);
            S3.iInsModeSeq = j3;
            
        case 2  % Radio mode change
            
            id = RadModeSeq(j3).id;
            [S3.RadModeDescr, S3.iRadModeDescr] = get_mode_descr(RadModeDescrList, id);
            S3.iRadModeSeq = j3;
            
        case 3  % Classification (selection/rejection) of (unclassified) rich data
            
            Classif = ClassifSeq(j3);
            
            % ASSERTIONS
            if (Classif.selectedRichBytes < 0) || (Classif.rejectedRichBytes < 0)
                error('Illegally configured rich data classification event: negative amount of bytes.')
            end
            if (Classif.rejectedRichBytes + Classif.selectedRichBytes) > S3.unclasRichBytes
                error('Illegally configured classification event: rejecting and selecting (sum) more data than available.')
                %Classif.rejectedRichBytes = S3.unclasRichBytes / 2;
                %Classif.selectedRichBytes = S3.unclasRichBytes / 2;
            end
            
            S3.unclasRichBytes = S3.unclasRichBytes - (Classif.selectedRichBytes + Classif.rejectedRichBytes);
            S3.queuedRichBytes = S3.unclasRichBytes + (Classif.selectedRichBytes);
            
            % NOTE: Event could lead to
            %   (1) on storage limit --> not on storage limit
            %   (2) ... --> reaching zero bytes (unclasRichBytes)
            % Could be important for setting flags.
            
        case 4  % Downlink change

            % ASSERTION. Remove?
            if DownlinkSeq(j3).bandwidthBps < 0
                error('Illegal negative downlink bandwidth.')
            end
            
            S3.downlinkBps = DownlinkSeq(j3).bandwidthBps;
            S3.iDownlinkSeq = j3;
    end
    iCommEvent = iCommEvent + 1;
    
    timeSec2 = timeSec3;
    S2 = S3;
end    % while

S1 = S3;

StorageState1 = [];
StorageState1.unclasRichBytes = S1.unclasRichBytes;
StorageState1.queuedRichBytes = S1.queuedRichBytes;
StorageState1.queuedSurvBytes = S1.queuedSurvBytes;

Clf.surv = SystemPrps.survClf;
Clf.rich = SystemPrps.richClf;

end



function [state1, timeSecP, StateArrays] = evolve_state_linearly(SystemPrps, state0, timeSec0, timeSec1, StateArrays)
% Evolve (integrate, extrapolate) the state of system over time period when there is no commanded event.
% This function defines very much of how the system evolves over time and is as such central.
% NOTE: This function defines which bounds (maximum storage, amount of data always positive) are followed.
% NOTE: This function defines how data is prioritized oat downlinking and destruction.
%
% ARGUMENTS
% =========
% state0, timeSec0   : Initial state at timeSec1.
% timeSec1           : Time to which the function will attempt to evolve the state.
%
% RETURN VALUES
% =============
% state1   : State at timeSec1. If timeSecP < timeSec1, then empty.
% timeSecP : The time until which evolution was successful. timeSecP <= timeSec1.
%
% BUG?: Numerical inexactness could lead to setting non-zero destroy* variables which is NOT enough to make usedStorageBps == 0.

% TODO: Add state arrays.
% TODO: New name? evolve_state, ~linear_step, ~step
%   NOTE: Not really a linear step since destroy* and excess downlink change abruptly (discontinuous).

% DEBUG
%timeSec0
%timeSec1

S0 = state0; clear state0
Sa = StateArrays; clear StateArrays

%==================================
% Initialize the state at timeSec0
%=================================
S0 = complement_state(S0, SystemPrps);


%===============================================
% Evolve non-constant state variables over time
%
% Use linear extrapolation to evolve those variables which values are determined by bounds.
% Check if they exceed those bounds.
%===============================================

% Reach maximum storage limit?
% NOTE: S1.usedStorageBytes is calculated from other variables (rather than integrated over time here; this avoids
% errors/inconsistencies with other variables).
S0.usedStorageBytes = min(S0.usedStorageBytes, SystemPrps.storageBytes);
[timeStorageMaxSec, ~, Sa.usedStorageBytes] = linear_extrapolate_limit(...
    timeSec0, S0.usedStorageBytes, timeSec1, S0.usedStorageBps/8, ...
    Sa.timeArraySec, Sa.usedStorageBytes, -Inf, SystemPrps.storageBytes);
S1 = [];

% Reach zero unclassified rich data?
%unclasRichBytes0     = S0.unclasRichBytes         % DEBUG
%changeUnclasRichBps0 = S0.changeUnclasRichBps     % DEBUG
S0.unclasRichBytes = max(S0.unclasRichBytes, 0);
[timeUnclasifRichMinSec, S1.unclasRichBytes, Sa.unclasRichBytes] = linear_extrapolate_limit(...
    timeSec0, S0.unclasRichBytes, timeSec1, S0.changeUnclasRichBps/8, ...
    Sa.timeArraySec, Sa.unclasRichBytes, 0, Inf);

% Reach zero queued rich data?
S0.queuedRichBytes = max(S0.queuedRichBytes, 0);
[timeQueuedRichMin, S1.queuedRichBytes, Sa.queuedRichBytes] = linear_extrapolate_limit(...
    timeSec0, S0.queuedRichBytes, timeSec1, S0.changeQueuedRichBps/8, ...
    Sa.timeArraySec, Sa.queuedRichBytes, 0, Inf);

% Reach zero survey data?
[timeQueuedSurvMinSec, S1.queuedSurvBytes, Sa.queuedSurvBytes] = linear_extrapolate_limit(...
    timeSec0, S0.queuedSurvBytes, timeSec1, S0.changeQueuedSurvBps/8, ...
    Sa.timeArraySec, Sa.queuedSurvBytes, 0, Inf);

timeSecP = min([timeSec1, timeStorageMaxSec, timeUnclasifRichMinSec, timeQueuedRichMin, timeQueuedSurvMinSec], [], 'omitnan');



%============================================================
% Handle variables which are constant over the time interval
%============================================================
CONSTANT_STATE_VARIABLES = {'downlinkBps', 'prodSurvBps', 'prodRichBps', 'iInsModeDescr', 'iRadModeDescr', 'powerWatt'};
%CONSTANT_ARRAY_FIELDS = {'iInsModeDescr', 'iRadModeDescr', 'prodSurvBps', 'prodRichBps', 'powerWatt'};
CONSTANT_ARRAY_FIELDS = CONSTANT_STATE_VARIABLES;
for i = 1:numel(CONSTANT_ARRAY_FIELDS)
    fn = CONSTANT_ARRAY_FIELDS{i};
    Sa.(fn) = fill_array_interval(timeSec0, timeSec1, S0.(fn), Sa.timeArraySec, Sa.(fn));
end
% Copy state variables that have not changed.
%CONSTANT_STATE_FIELDS = {'downlinkBps', 'prodSurvBps', 'prodRichBps', 'InsModeDescr', 'RadModeDescr', 'iInsModeDescr', 'iRadModeDescr'};
CONSTANT_STATE_FIELDS = {CONSTANT_STATE_VARIABLES{:}, 'InsModeDescr', 'RadModeDescr'};
for i = 1:numel(CONSTANT_STATE_FIELDS)
    fn = CONSTANT_STATE_FIELDS{i};
    S1.(fn) = S0.(fn);
end

% Should not be necessary. Just for safety.
if timeSecP < timeSec1
    S1 = [];
end



state1 = S1;
StateArrays = Sa;
end



function State = complement_state(State, SystemPrps)
% Using existing state variables to derive additional state variables, valid at the same instant:
%   (1) downlink distributed on data types
%   (2) destroy*
%   (3) usedStorageBps
%   (4) total survey and rich data production rates
%   (5) total power

S = State;

% Compile total production of survey & rich data respectively.
S.prodSurvBps = S.InsModeDescr.prodSurvBps + S.RadModeDescr.prodSurvBps;
S.prodRichBps = S.InsModeDescr.prodRichBps + S.RadModeDescr.prodRichBps;
S.powerWatt   = S.InsModeDescr.powerWatt   + S.RadModeDescr.powerWatt;


% Choose how downlink is distributed over different types of stored data.
[S.downlinkSurvBps, S.downlinkRichBps, S.downlinkExceBps] = distribute_value(...
    S.downlinkBps, ...
    S.queuedSurvBytes <= 0, S.prodSurvBps, ...
    S.queuedRichBytes <= 0, S.prodRichBps);

S.usedStorageBytes = ...
    SystemPrps.survClf * (S.queuedSurvBytes) + ...
    SystemPrps.richClf * (S.queuedRichBytes + S.unclasRichBytes);

S.destroyQueuedSurvBps = 0;
S.destroyUnclasRichBps = 0;
S.destroyQueuedRichBps = 0;   % Longer variable name just to keep it analogous with other destroy* variables.

S = calc_used_storage_and_change_bps(S, SystemPrps);

% If used storage about to exceed storage limit, then set destroy* variables to cancel out.
if (S.usedStorageBytes >= SystemPrps.storageBytes) && (S.usedStorageBps > 0)
    % CASE: Used storage about to exceed storage limit.

    % Set destroy* values depending on usedStorageBps.
    S.destroyStorageBps = S.usedStorageBps;
    
    % Choose how data destruction is distributed over different types of stored data.
    [destroyUnclasRichStorageBps, destroyQueuedRichStorageBps, destroySurvStorageBps] = distribute_value(...
       S.destroyStorageBps, ...
        S.unclasRichBytes <= 0, S.prodRichBps*SystemPrps.richClf, ...
        S.queuedRichBytes <= 0, 0);
    
    S.destroyUnclasRichBps = destroyUnclasRichStorageBps / SystemPrps.richClf;
    S.destroyQueuedRichBps = destroyQueuedRichStorageBps / SystemPrps.richClf;
    S.destroyQueuedSurvBps = destroySurvStorageBps       / SystemPrps.survClf;
    
    % ASSERTION
    S = calc_used_storage_and_change_bps(S, SystemPrps);
    if S.usedStorageBps ~= 0
        error('usedStorageBps ~= 0 after setting destroy* variables.')
    end
end

State = S;
end



function State = calc_used_storage_and_change_bps(State, SystemPrps)
S = State;

S.changeQueuedSurvBps = S.prodSurvBps - S.downlinkSurvBps - S.destroyQueuedSurvBps;
S.changeUnclasRichBps = S.prodRichBps                     - S.destroyUnclasRichBps;
S.changeQueuedRichBps =               - S.downlinkRichBps - S.destroyQueuedRichBps;

S.usedStorageBps = ...
    SystemPrps.survClf * (S.changeQueuedSurvBps) + ...
    SystemPrps.richClf * (S.changeUnclasRichBps + S.changeQueuedRichBps);

State= S;
end


    
% Distribute (split, spread) a finite positive value over three variables.
%
% Prioritize val1, then val2, then val3.
% If corresponding conditions maxCond1, maxCond2 are true, then the corresponding variable can only be set to the
% corresponding positive maximum value.
%
% val1+val2+val3 == value
function [val1, val2, val3] = distribute_value(value, maxCond1, max1, maxCond2, max2)
    if maxCond1
        val1 = min(max1, value);
    else
        val1 = value;
    end
    value = value - val1;
    if maxCond2
        val2 = min(max2, value);
    else
        val2 = value;
    end
    value = value - val2;
    val3 = value;
end



% Set parts of array to constant value.
% Should be identical to linear_extrapolate_limit with dydx==0, yMin==-Inf, yMax==Inf.
function [yArray] = fill_array_interval(x0, x1, y0, xArray, yArray)
i = (x0 <= xArray) & (xArray <= x1);
yArray(i) = y0;
end



function [xp, y1, yArray] = linear_extrapolate_limit(x0, y0, x1, dydx, xArray, yArray, yMin, yMax)
% Assume function y=y(x). Linearly extrapolate y values from x0,y0 to x1,y1 but keep track of when exceeding bounds.
%
% ARGUMENTS
% =========
% x0, y0     : Starting point for extrapolation
% dydx       : Constant derivative used for the extrapolation.
% x1              : Maximum end point for extrapolation.
% yMin, yMax      : Interval of y values which the extrapolation x0,y0 to xp,yp must stay within.
% xArray, yArray  : Arrays with x and (potentially) y values.
%
%
% RETURN VALUES
% =============
% xp     : Point where extrapolation reached min or max limit. If y0 is outside the yMin-yMax interval, then xp==x0,
%          yp==y0, i.e. yp is OUTSIDE the yMin-yMax interval for this special case.
% y1     : y value for x=x1. Can be outside bounds.
% yArray : Argument yArray, but with the y values for the range x0-x1 (in xArray) set to yArray(i)==y(xArray(i)).
%
%

% TODO-DECISION: Inconsistent to NOT treat y0 out-of-bounds as special case, with special flag?
%
% PROPOSAL: Add flag for whether y0 is outside bounds.
% PROPOSAL: Add assertions.

% Set return values which apply regardless of min-max boundaries.
i = (x0 <= xArray) & (xArray <= x1);
yArray(i) = y0 + (xArray(i) - x0) * dydx;
y1        = y0 + (x1        - x0) * dydx;

%============
% Set xp, yp
%============
% Check for y0 being out of bounds.
[~, reachedMin, reachedMax] = outside_boundaries(y0, yMin, yMax);
if (reachedMin || reachedMax)
    xp = x0;
    return
end
% Check for y0 < y <= y1 being out of bounds (knowing that we are working with a linear function y=y(x)).
[yp, reachedMin, reachedMax] = outside_boundaries(y1, yMin, yMax);
if (reachedMin || reachedMax)
    % NOTE: Does not work for dydx==0, but that should never happen here, since
    % dydx==0   ==>   y0==y1   ==>   out-of-boundaries for y1, if ever.
    xp = x0 + (yp - y0) / dydx;
else
    xp = x1;
end
end



% yp : If inside yMin-yMax, then equal to y, otherwise that value yMin, yMax which is closest/exceeded.
function [yp, belowMin, aboveMax] = outside_boundaries(y, yMin, yMax)
belowMin = false;
aboveMax = false;

if y < yMin
    belowMin = true;
    yp       = yMin;
elseif y > yMax
    aboveMax = true;
    yp       = yMax;
else
    yp = y;
end
end



function [ModeDescr, iModeDescr] = get_mode_descr(ModeDescrList, id)
iModeDescr = find(strcmp(id, {ModeDescrList.id}));
if numel(iModeDescr) ~= 1
    error('Can not find exactly one mode description for a given mode ID.')
end
ModeDescr = ModeDescrList(iModeDescr);
end



% Return factor describing how much onboard storage space is used for a file of a given size.
% clusterLossFactor = <bytes allocated on storage> / <bytes of data stored in file>   >=   1
function clusterLossFactor = cluster_loss_factor(fileDataSize, clusterSize)
nClusters         = fileDataSize / clusterSize;
fileStorageUse    = ceil(nClusters) * clusterSize;
clusterLossFactor = fileStorageUse / fileDataSize;
end



% In an array y, find the last value that is smaller than or equal to y0.
function [i] = prev_lower_equal_value(y, y0)
i  = find(y <= y0, 1, 'last');
%yl = y(i);
end



% Generic utility function.
% Merge and sort list of 1D vectors, and return the indices needed to find the corresponding source elements.
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% varargin : Sequence of numeric 1D arrays. The arrays do not need to be sorted.
% y        : Numeric 1D array that is the sorted concatenation of the varargin arrays.
% i,j      : Numeric 1D arrays such that y(n) == varargin{i(n)}(j(n)), for all n.
function [y, i, j] = merge_seq(varargin)
yList = {};
iList = {};
jList = {};
for k = 1:numel(varargin)
    temp     = varargin{k};
    
    % ASSERTION
    if ~isnumeric(temp)
        error('Argument is not numeric')
    end
    
    yList{k} = temp(:)';   % Force to be row vector.
    
    iList{k} = ones(1, numel(yList{k})) * k;   % Row vector [k, k, ..., k].
    jList{k} = 1:numel(yList{k});              % Row vector [1, 2, ...   ].
end

yUnsorted = [yList{:}];   % Merge row vectors into one.
iUnsorted = [iList{:}];   % Merge row vectors into one.
jUnsorted = [jList{:}];   % Merge row vectors into one.

[y, kSort] = sort(yUnsorted);
i = iUnsorted(kSort);
j = jUnsorted(kSort);
end
