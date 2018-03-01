function [StorageState1, StateArrays, Clf] = engine(Constants, StorageState0, CommEvents, timeSec0, timeSec1, timeArraySec)
% Engine for JUICE TM/power budget simulation.
%
%
%   IMPORTANT NOTE: Classification events have been disabled. Rich data production is directed to
%                   queued rich data. CommEvents.ClassifSeq has been disabled and must be set to [].
%
%
% PROBLEM INTENDED TO BE SOLVED
% =============================
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
% DESIGN GOALS
% ============
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
% Clf                   : Struct with cluster loss factors.
%   .surv
%   .rich
%
%
% VARIABLE NAMING CONVENTIONS, DEFINITIONS OF TERMS
% ===============================================
% bps = Bits per second (not BYTES per second)
% Int = Integrated over time
% Rad = Radio (mode)
% Ins = In situ (mode)
% Surv = Survey (data)
% Rich = Rich (data)
% Seq = Sequence
% CLF = Cluster Loss Factor = The ratio between allocated storage and data in file (>=1).
%
%
% NOTES
% =====
% ** Not running any instrument mode is modelled by defining and "running" a made-up "No-instrument-mode instrument
% mode" which produces zero data.
% ** Not having any downlink is modelled by setting Downlink(i).bandwidthBps == 0.
% ** Only includes the power consumption from the in situ and radio modes, not from "powerDPU = 1912; powerXtra = 1100;"
% as in "ju_rpwi_tm_power.m".
% ** Classification events are handled if-and-only-if they occur at times timeSec0 <= timeSec < timeSec1.
%
%
% IMPLEMENTATION NOTES
% ====================
% ** Classification events instantly influence the storage state (which is returned) and must therefore be handled
% properly at the ends of the time interval to avoid executing the same classification events twice or never when
% chaining calls to this function (indirectly through the main function "engine"), one after another to evolve a longer
% time interval. The implementation MUST therefore EITHER execute classification events at (1) timeSec0, OR (2)
% timeSec1, but not at both, nor neither. Note also behaviour for the special case of zero-length time interval.
% ** Amount of stored data is specified in bytes since MB, GB etc would be slightly ambiguous.
% ** Only permitting two file sizes, for survey and rich data separately, to keep the model simple. Having more file
% sizes (e.g. one per mode) complicates the representation of stored data since the relationship between (1) amount of
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
%
% PROPOSAL: evolve --> integrate
%
% PROPOSAL: Separate name for "State" (Storage State+mode+downlink), to distinguish from "Storage State".
%   PROPOSAL: SystemState
%       CON: Will have identical shortening, "SS"
%           PROPOSAL: "SyS" vs "StS".
%   PROPOSAL: InstrumentState.
%       CON: Could be interpreted as only instrument modes (no downlink, no storage state; storage could be seen as
%       separate from instrument)
%   PROPOSAL: Flow State
%   PROPOSAL: StorageState-->DataState
%       CON: Does not imply anything with the actual drive on the s/c. Too generic.
%
% PROPOSAL: Make sure to have variables for total (integrated) amount of data generated for any data type, regardless of whether some of it has been destroyed or not.
%   NOTE: Basically possible with the old version of engine, and with step_func.
%   PRO: Useful for checking realism of commanded modes.
% 
% PROPOSAL: Have function chose the timestamps.
%   PROPOSAL: Choose only (1) all discrete events, and (2) beginning and end.
%   PROBLEM: How handle (true) step functions?


% Internal variable naming conventions
% ------------------------------------
% S      = State (not storage state)
% CS     = Complemented state (a state + derived variables)
% 0, 1   = Initial and final state/time for call to function that evolves state.
% 0b, 1b = Like 0, 1 but for a subset of the evolution.
% CLF    = Cluster Loss Factor (>=1)



% ASSERTION
if ~isempty(CommEvents.ClassifSeq)
    error('engine:Assertion', 'CommEvents.ClassifSeq is not empty. NOTE: Classification events have been disabled. Rich data production is directed to queued rich data.')
end
CommEvents.ClassifSeq = struct('timeSec',  {}, 'selectedRichBytes', {}, 'rejectedRichBytes', {});



% Initialize StateArrays.
% Initialize fields with NaN to be able to detect unset values before exiting.
EMPTY_ARRAY_FIELDS = {'unclasRichBytes', 'queuedRichBytes', 'queuedSurvBytes', 'iInsModeDescr', 'iRadModeDescr', ...
    'prodSurvBps', 'prodRichBps', 'powerWatt', 'usedStorageBytes', 'downlinkBps', 'downlinkSurvBps', 'downlinkRichBps', 'downlinkExceBps'};  % State array fields to initialize.
tempEmptyArray = zeros(size(timeArraySec)) * NaN;
StateArrays = [];
StateArrays.timeArraySec = timeArraySec;
for i = 1:numel(EMPTY_ARRAY_FIELDS)
    StateArrays.(EMPTY_ARRAY_FIELDS{i}) = tempEmptyArray;
end

% Derive CLFs.
Constants.SystemPrps.survClf = cluster_loss_factor(Constants.SystemPrps.survFileSizeBytes, Constants.SystemPrps.clusterSizeBytes);
Constants.SystemPrps.richClf = cluster_loss_factor(Constants.SystemPrps.richFileSizeBytes, Constants.SystemPrps.clusterSizeBytes);



%====================
% Initialize state 0
%====================
State0 = StorageState0;   % Incomplete as is. Fields will be added to it.
% Find initial mode and downlink, by finding the corresponding most recent timeSec <= timeSec0 events.
State0.iInsModeSeq  = prev_lower_equal_value([CommEvents.InsModeSeq.beginSec],  timeSec0);
State0.iRadModeSeq  = prev_lower_equal_value([CommEvents.RadModeSeq.beginSec],  timeSec0);
State0.iDownlinkSeq = prev_lower_equal_value([CommEvents.DownlinkSeq.beginSec], timeSec0);
% ASSERTIONS
if isempty(State0.iInsModeSeq)
    error('Can not find in-situ mode covering beginning of time interval.')
elseif isempty(State0.iRadModeSeq)
    error('Can not find radio mode covering beginning of time interval.')
elseif isempty(State0.iDownlinkSeq)
    error('Can not find downlink covering beginning of time interval.')
end



[State1, StateArrays] = evolve_state(Constants, CommEvents, State0, StateArrays, timeSec0, timeSec1);



% Construct additional return values.
StorageState1 = [];
StorageState1.unclasRichBytes = State1.unclasRichBytes;
StorageState1.queuedRichBytes = State1.queuedRichBytes;
StorageState1.queuedSurvBytes = State1.queuedSurvBytes;
%
Clf.surv = Constants.SystemPrps.survClf;
Clf.rich = Constants.SystemPrps.richClf;

end



function [State1, StateArrays] = evolve_state(Constants, CommEvents, State0, StateArrays, timeSec0, timeSec1)
% Evolve the state over arbitrary time period.
%
% This function 
% (1) calls evolve_state_wo_CE.
% (2) handles commanded events
%
% NOTE: Classification events are handled if-and-only-if they occur at times timeSec0 <= timeSec < timeSec1.
% NOTE: The most recent non-classification events at timeSec<timeSec0 are handled outside of this function in order to
% initialize the submitted state.
% 

assert_state(State0);
%fprintf('evolve_state: %g -- %g\n', timeSec0, timeSec1)   % DEBUG

% Create one single, time-ordered list of commanded events.
[commEventTimesSec, iEventType, jEventTypeSeq] = merge_seq(...
    [CommEvents.InsModeSeq.beginSec], ...
    [CommEvents.RadModeSeq.beginSec], ...
    [CommEvents.ClassifSeq.timeSec], ...
    [CommEvents.DownlinkSeq.beginSec]);

% Find first commanded event.
iCommEvent = find(commEventTimesSec >= timeSec0, true, 'first');


timeSec0b = timeSec0;
State0b   = State0;

noMoreCommEvent = false;
while true
    %timeSec0b
    %State0b   % DEBUG

    %================================================================
    % Set timeSec1b
    % -------------
    % Set timeSec1b to be the value that comes first of
    % (1) end of interval,
    % (2) the next commanded event.
    %
    % NOTE: May have timeSec0b == timeSec1b if
    % (1) there is a commanded event at timeSec0, or
    % (2) there are multiple commanded events at the same timestamp.
    %================================================================
    if (iCommEvent > numel(commEventTimesSec)) || (commEventTimesSec(iCommEvent) > timeSec1)
        % CASE: (1) There is no more (caller-specified) commanded event, or
        %       (2) the next commanded event takes place after timeSec1.
        noMoreCommEvent = true;
        timeSec1b = timeSec1;
    else
        % CASE: There is at least one relevant, future commanded event.
        % ==> Evolve until it.
        timeSec1b = commEventTimesSec(iCommEvent);
    end

    % Evolve state 0b-->1b
    [State1b, StateArrays] = evolve_state_wo_CE(Constants, CommEvents, State0b, StateArrays, timeSec0b, timeSec1b);
    
    %timeSec1b
    %State1b   % DEBUG

    if (timeSec1b == timeSec1) && noMoreCommEvent
        State0b = State1b;   % Value that will be read outside of loop, and which must be initialized if zero iterations.
        break   % EXIT LOOP
    end

    % ASSERTION
    if timeSec1b > timeSec1
        error('engine:Assertion', 'timeSec1b > timeSec1')
    end

    %=============================================
    % Handle ONE commanded event (modify State1b)
    %=============================================
    jSeq = jEventTypeSeq(iCommEvent);
    switch iEventType(iCommEvent)
        
        case 1  % In situ mode change
            
            State1b.iInsModeSeq = jSeq;
            
        case 2  % Radio mode change
            
            State1b.iRadModeSeq = jSeq;
            
        case 3  % Classification (selection/rejection) of (unclassified) rich data
            
            Classif = CommEvents.ClassifSeq(jSeq);
            
            % ASSERTIONS
            % NOTE: Assertion for checking the configuration of the classification event, NOT the final state.
            % The final state is checked via assert_state (not max storage limit).
            if (Classif.selectedRichBytes < 0) || (Classif.rejectedRichBytes < 0)
                error('Illegally configured rich data classification event: negative amount of bytes.')
            end
            if (Classif.rejectedRichBytes + Classif.selectedRichBytes) > State1b.unclasRichBytes
                error('Illegally configured classification event: rejecting and selecting (sum) more data than available.')
                %Classif.rejectedRichBytes = S3.unclasRichBytes / 2;
                %Classif.selectedRichBytes = S3.unclasRichBytes / 2;
            end

            State1b.unclasRichBytes = State1b.unclasRichBytes - (Classif.selectedRichBytes + Classif.rejectedRichBytes);
            State1b.queuedRichBytes = State1b.queuedRichBytes + (Classif.selectedRichBytes);

            % NOTE: Event could lead to
            %   (1) on storage limit --> not on storage limit
            %   (2) ... --> reaching zero bytes (unclasRichBytes)
            % Could be important for setting flags.
            
        case 4  % Downlink change
    
            State1b.iDownlinkSeq = jSeq;
    end

    iCommEvent = iCommEvent + 1;

    timeSec0b = timeSec1b;
    State0b   = State1b;
    assert_state(State0b);
    
end    % while
% CASE: State0b contains the latest state. NOTE: Must also work if no loop iteration!

State1 = State0b;
assert_state(State1);
assert_state_arrays(StateArrays);

end



function [State1, StateArrays] = evolve_state_wo_CE(Constants, CommEvents, State0, StateArrays, timeSec0, timeSec1)
% Evolve the state over time period where
% (1) there is no commanded event.
%
% In practice, this function is a wrapper around evolve_state_wo_CE_discont. This function will chain calls to that
% function in order to ALWAYS COMPLETE the entire specified time period.
%
% NOTE: Permits zero-length time period.
%


assert_state(State0)
%fprintf('evolve_state_wo_CE: %g -- %g\n', timeSec0, timeSec1)   % DEBUG

timeSec0b = timeSec0;
State0b   = State0;
while timeSec0b < timeSec1      % timeSec0 == timeSec1 ==> No iteration, avoid calling evolve_state_wo_CE_discont.
    % TRY to evolve system for 0b-->1.
    [State1b, timeSec1b, StateArrays] = evolve_state_wo_CE_discont(Constants, CommEvents, ...
        State0b, StateArrays, timeSec0b, timeSec1);
    % SUCCEEDED in evolving system for 0b-->1b.
    
    if timeSec1b < timeSec1
        % CASE: Did not evolve until 1, only part of the way.
        % ==> Try again for the remaining time period.
        timeSec0b = timeSec1b;
        State0b   = State1b;
    else
        % CASE: Succeeded in evolving until 1.
        % ==> Quit
        State0b = State1b;
        break
    end
end
% CASE: State0b is the latest state. NOTE: Must also work with zero iterations.

State1 = State0b;
end



function [StateP, timeSecP, StateArrays] = evolve_state_wo_CE_discont(Constants, CommEvents, State0, StateArrays, timeSec0, timeSec1)
% Evolve (integrate, extrapolate) the state of system over time period
% (1) without commanded events (CE), and
% (2) without discontinuities (discont), i.e. system evolves analytically (linearly) without sudden jumps due to
% reaching limits.
%
% The function stops early when it can not satisfy (2) (interrupts when reaching zero bytes in any data category; when
% reaching storage max limit).
%
% This function defines very much of how the system evolves over time and is as such central.
%
% NOTE: This function defines which bounds (maximum storage, amount of data always positive) are followed.
% NOTE: This function defines how data is prioritized at downlinking and destruction (producing too much data).
%
%
% ARGUMENTS
% =========
% State0    : Initial state at timeSec0.
% timeSec0
% timeSec1  : Time to which the function will TRY to evolve the state.
%
%
% RETURN VALUES
% =============
% StateP   : State at timeSecP.
% timeSecP : The time until which the evolution was successful. timeSecP <= timeSec1.

%fprintf('evolve_state_wo_CE_discont: %.18g -- %.18g\n', timeSec0, timeSec1)   % DEBUG

assert_state(State0);

S0 = State0;      clear State0
Sa = StateArrays; clear StateArrays

CS0 = complement_state(S0, Constants, CommEvents);   % CS = Complemented State.



timeSecP = timeSec1;
SP = [];
nbrOfTries = 0;
while true
    fprintf('Try 4 * linear_extrapolate_limit: %.18g -- %.18g\n', timeSec0, timeSecP);    % DEBUG

    % ASSERTION
    % Prevent algorithm from getting stuck in a loop. Wrapper function should handle case of evolving over a zero-length
    % time interval so this function should never have to.
    %if (timeSec0 == timeSecP)
    %    error('evolve_state_wo_CE_discont:Assertion', 'linear_extrapolate_limit:Assertion', 'timeSec0 <= timeSecP')
    %end
    % NOTE: linear_extrapolate_limit requires non-zero time period (assertion). -- NOT?!

    % Evolve ALLOCATED STORAGE.
    % NOTE: This purpose is
    % * NOT to update any actual state variable
    % * to check if allocated storage is within limits
    % * to update a corresponding state array
    [timeStorageMaxSec, ~, Sa.usedStorageBytes] = TM_power_budget.engine_utils.linear_extrapolate_limit(...
        timeSec0, timeSecP, CS0.usedStorageBytes, CS0.usedStorageBps/8, ...
        Sa.timeArraySec, Sa.usedStorageBytes, -Inf, Constants.SystemPrps.storageBytes);

    % Evolve UNCLASSIFIED RICH data
    [timeUnclasifRichMinSec, SP.unclasRichBytes, Sa.unclasRichBytes] = TM_power_budget.engine_utils.linear_extrapolate_limit(...
        timeSec0, timeSecP, CS0.unclasRichBytes, CS0.changeUnclasRichBps/8, ...
        Sa.timeArraySec, Sa.unclasRichBytes, 0, Inf);

    % Evolve QUEUED RICH data.
    [timeQueuedRichMinSec, SP.queuedRichBytes, Sa.queuedRichBytes] = TM_power_budget.engine_utils.linear_extrapolate_limit(...
        timeSec0, timeSecP, CS0.queuedRichBytes, CS0.changeQueuedRichBps/8, ...
        Sa.timeArraySec, Sa.queuedRichBytes, 0, Inf);
    
    % Evolve (QUEUED) SURVEY data.
    [timeQueuedSurvMinSec, SP.queuedSurvBytes, Sa.queuedSurvBytes] = TM_power_budget.engine_utils.linear_extrapolate_limit(...
        timeSec0, timeSecP, CS0.queuedSurvBytes, CS0.changeQueuedSurvBps/8, ...
        Sa.timeArraySec, Sa.queuedSurvBytes, 0, Inf);

    nbrOfTries = nbrOfTries + 1;
    timeSecP = min(...
        [timeStorageMaxSec, timeUnclasifRichMinSec, timeQueuedRichMinSec, timeQueuedSurvMinSec], ...
        [], 'omitnan');   % The time until which the evolution was successful.
    if (timeSecP == timeSec1) || (nbrOfTries >=2)
        break
    end
end

SP.iInsModeSeq  = S0.iInsModeSeq;
SP.iRadModeSeq  = S0.iRadModeSeq;
SP.iDownlinkSeq = S0.iDownlinkSeq;



%======================================================
% Create state arrays for variables which are constant
%======================================================
CONSTANT_STATE_ARRAYS = {'downlinkBps', 'downlinkSurvBps', 'downlinkRichBps', 'downlinkExceBps', 'prodSurvBps', 'prodRichBps', 'iInsModeDescr', 'iRadModeDescr', ...
    'powerWatt'};
for i = 1:numel(CONSTANT_STATE_ARRAYS)
    fn = CONSTANT_STATE_ARRAYS{i};
    Sa.(fn) = fill_array_interval(timeSec0, timeSec1, CS0.(fn), Sa.timeArraySec, Sa.(fn));
end



StateP = SP;
StateArrays = Sa;
assert_state(StateP);
end



function ComplementedState = complement_state(State, Constants, CommEvents)
% Using existing state variables to derive additional state variables, valid at the same instant:
%   (1) downlink distributed on data types
%   (2) destroyQueuedSurvBps
%       destroyUnclasRichBps
%       destroyQueuedRichBps
%   (3) usedStorageBps
%   (4) total survey and rich data production rates
%   (5) total power

assert_state(State);

CS = State; clear State   % CS = Complemented State
SystemPrps = Constants.SystemPrps;



% Look up values in CommEvents.
[CS.iInsModeDescr, InsModeDescr] = get_mode_descr(CommEvents.InsModeSeq, Constants.InsModeDescrList, CS.iInsModeSeq);
[CS.iRadModeDescr, RadModeDescr] = get_mode_descr(CommEvents.RadModeSeq, Constants.RadModeDescrList, CS.iRadModeSeq);
CS.downlinkBps = CommEvents.DownlinkSeq(CS.iDownlinkSeq).bandwidthBps;

% Compile total production of survey & rich data respectively + power consumption.
CS.prodSurvBps = InsModeDescr.prodSurvBps + RadModeDescr.prodSurvBps;
CS.prodRichBps = InsModeDescr.prodRichBps + RadModeDescr.prodRichBps;
CS.powerWatt   = InsModeDescr.powerWatt   + RadModeDescr.powerWatt;

% Choose how downlink is distributed over different types of stored data.
if 0
    [CS.downlinkSurvBps, CS.downlinkRichBps, CS.downlinkExceBps] = distribute_value(...
        CS.downlinkBps, ...
        CS.queuedSurvBytes <= 0, CS.prodSurvBps, ...
        CS.queuedRichBytes <= 0, 0);    % NOTE: No QUEUED rich data is produced. Can therefore not be downlinked if zero.
else
    [CS.downlinkSurvBps, CS.downlinkRichBps, CS.downlinkExceBps] = distribute_value(...
        CS.downlinkBps, ...
        CS.queuedSurvBytes <= 0, CS.prodSurvBps, ...
        CS.queuedRichBytes <= 0, CS.prodRichBps);    % NOTE: QUEUED rich data is produced.
end

CS.usedStorageBytes = ...
    SystemPrps.survClf * (CS.queuedSurvBytes) + ...
    SystemPrps.richClf * (CS.queuedRichBytes + CS.unclasRichBytes);
CS.usedStorageBytes = min(CS.usedStorageBytes, SystemPrps.storageBytes);

CS.destroyQueuedSurvBps = 0;
CS.destroyUnclasRichBps = 0;
CS.destroyQueuedRichBps = 0;   % Longer variable name just to keep it analogous with other destroy* variables.

CS = set_change_fields(CS);

CS.usedStorageBps = ...
   SystemPrps.survClf * (CS.changeQueuedSurvBps) + ...
   SystemPrps.richClf * (CS.changeUnclasRichBps + CS.changeQueuedRichBps);

% If used storage about to exceed storage limit, then set destroy* variables to cancel out.
if (CS.usedStorageBytes >= SystemPrps.storageBytes) && (CS.usedStorageBps > 0)
    % CASE: Used storage about to exceed storage limit.

    % Set destroy* values depending on usedStorageBps.
    CS.destroyStorageBps = CS.usedStorageBps;

    % Choose how data destruction is distributed over different types of stored data.
    if 0
        [destroyUnclasRichStorageBps, destroyQueuedRichStorageBps, destroySurvStorageBps] = distribute_value(...
            CS.destroyStorageBps, ...
            CS.unclasRichBytes <= 0, CS.prodRichBps*SystemPrps.richClf, ...
            CS.queuedRichBytes <= 0, 0);
    else
        [destroyUnclasRichStorageBps, destroyQueuedRichStorageBps, destroySurvStorageBps] = distribute_value(...
            CS.destroyStorageBps, ...
            CS.unclasRichBytes <= 0, 0, ...
            CS.queuedRichBytes <= 0, CS.prodRichBps*SystemPrps.richClf);
    end
    
    CS.destroyUnclasRichBps = destroyUnclasRichStorageBps / SystemPrps.richClf;
    CS.destroyQueuedRichBps = destroyQueuedRichStorageBps / SystemPrps.richClf;
    CS.destroyQueuedSurvBps = destroySurvStorageBps       / SystemPrps.survClf;

    % POTENTIAL BUG? Calculated values may be slightly wrong due to numerical inaccuracies?!! Exact values are required
    % (more or less for sensible assertions) for evolution of state?!
    CS = set_change_fields(CS)
    
    CS.usedStorageBps = 0;
    
    % ASSERTION. Can be triggered due to numerical inaccuracy.
    %if CS.usedStorageBps ~= 0
    %    error('usedStorageBps ~= 0 after setting destroy* variables.')
    %end
end

ComplementedState = CS;
end



function CS = set_change_fields(CS)
    CS.changeQueuedSurvBps = CS.prodSurvBps - CS.downlinkSurvBps - CS.destroyQueuedSurvBps;
    if 0
        CS.changeUnclasRichBps = CS.prodRichBps                      - CS.destroyUnclasRichBps;
        CS.changeQueuedRichBps =                - CS.downlinkRichBps - CS.destroyQueuedRichBps;
    else
        CS.changeUnclasRichBps =                                       CS.destroyUnclasRichBps;
        CS.changeQueuedRichBps = CS.prodRichBps - CS.downlinkRichBps - CS.destroyQueuedRichBps;
    end
end



% Distribute (split, spread) a finite positive value over three variables.
%
% Prioritize val1, then val2, then val3.
% If corresponding conditions maxCond1, maxCond2 are true, then the corresponding variable can not be set to anything
% higher than corresponding positive maximum value.
%
% val1+val2+val3 == value
%
% NOTE: If one needs to distribute over two variables instead of three, just set maxCond1=true and max1=0, and ignore val1.
%
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
%
% ARGUMENTS
% =========
% x0, x1    : Scalars. Min & max of specified interval of xArray elements (not indices). Indirectly specifies a set of
%             indicies in yArray (and xArray).
% xArray    : 1D array.
% yArray    : 1D array of same size as xArray.
function [yArray] = fill_array_interval(x0, x1, y0, xArray, yArray)
i = (x0 <= xArray) & (xArray <= x1);
yArray(i) = y0;
end



function [iModeDescr, ModeDescr] = get_mode_descr(ModeSeq, ModeDescrList, iModeSeq)
% ASSUMES: ModeSeq.id exists.

id = ModeSeq(iModeSeq).id;
iModeDescr = find(strcmp(id, {ModeDescrList.id}));
if numel(iModeDescr) ~= 1
    error('Can not find exactly one mode description for a given mode ID="%s".', id)
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



function assert_state(State)
% PROPOSAL: Check max storage limit.
%   CON: Requires more arguments: CLFs, maxStorageBytes
%   CON: That is a derived quantity. ==> Too "model dependent". Not a true assertion. Should be checked when/where the
%        quantity is derived, maybe.
%   CON: ~Too complicated (for an assertion).

STATE_FIELDS = {...
    'iRadModeSeq', 'iInsModeSeq', 'iDownlinkSeq', ...
    'unclasRichBytes', 'queuedRichBytes', 'queuedSurvBytes'};

if ~isempty(setxor(fieldnames(State), STATE_FIELDS))
    error('assert_state:Assertion', 'Illegal set of fields.')
end
if (State.unclasRichBytes < 0) || (State.queuedRichBytes < 0) || (State.queuedSurvBytes < 0)
    error('assert_state:Assertion', 'Data type has less than zero bytes stored.')
end
end



function assert_state_arrays(StateArrays)
fieldNamesList = fieldnames(StateArrays);

for i = 1:numel(fieldNamesList)
    fn = fieldNamesList{i};
    if any(isnan(StateArrays.(fn)))
        error('State arrays field "%s" contains at least one NaN.', fn)
    end
end
end
