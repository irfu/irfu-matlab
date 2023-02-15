%
% Given a list of input intervals on the real number line, a <= x <= b,
% (1) find all "output intervals". Output intervals are uninterrupted by the
%     beginning or end of an input interval.
%     NOTE: Zero-length intervals (a<=x<=b, a=b) contribute with two such
%     interruptions (i.e. at "a-" and "b+").
%           MIGHT BE INCORRECT. ZLIs might only contribute with one
%           interruption.
%     They can
%     thus be seen as "infinitesimally short" rather than zero-length.
%           
%     NOTE: Adjacent intervals (a<=x<=b and b<=x<=c) only contribute with one
%     interruption where the intervals touch (i.e. at b).
% (2) for every output interval, find the set of input intervals that overlaps
%     with it.
%     NOTE: For output intervals, the boundaries may or may not be included. It
%     is in principle unambiguous, but in practice one has to investigate by
%     looking at intervals before and after.
%
%
% USAGE NOTE
% ==========
% To exclude that input interval boundaries count as overlapping, remove output
% sets for which h2Array-h1Array == 0.
% To exclude output intervals which are empty on input intervals, remove output
% sets for which nArray == 0.
%
%
% EDGE CASES (INPUT INTERVALS)
% ============================
%  ZLI =     Zero-Length Interval
% NZLI = Non-Zero Length Interval
% --
% Touching NZLIs count as overlapping
%   Ex: 1--2, 2--3
% ZLIs overlap with other ZLIs
%   Ex: 2--2, 2--2
% ZLIs overlap when they touch NZLIs
%   Ex: 2--2, 2--3
%   Ex: 1--2, 2--2, 2--3
%
%
% ARGUMENTS
% =========
% f1Array, 
% f2Array
%       1D numeric arrays, finite. Same length. intArray1(iInterval) and
%       intArray2(iInterval) define beginning and end of interval.
%       intArray1 =< intArray2. Intervals do not need to be sorted.
%
%
% RETURN VALUE
% ============
% setsCa
%       1D cell array (CA) of numeric 1D arrays. Indices into f1Array/f2Array.
%       {iInputInterval} = 1D numeric array. Sorted from low to high (adjacent)
%       intervals.
% nArray
%       nArray(iOutputInterval) = Number of input intervals in output interval.
% h1Array,
% h2Array
%       h1Array(iOutputInterval) < x < h2Array(iOutputInterval) define an output
%       interval.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-04-25.
%
function [setsCa, nArray, h1Array, h2Array] = find_interval_overlaps(f1Array, f2Array)
    % PROPOSAL: Make more consistent. Now supports two different kinds of input intervals: x=a, and a<=x<=b.
    %   Seem inconsistent somehow.
    %   PROPOSAL: Forbid zero-length input & output intervals. See input and output as having not overlapping boundaries (e.g. [a,b) ).
    %       CON: Can still be useful to know whether input intervals are touching (zero-length output intervals).
    %       PROPOSAL: Have other function with superset of functionality:
    %           Input & output intervals have flags for whether lower and upper boundary is included (separately for
    %           each interval, and for each boundary, i.e. 2N flags)
    %           [setsCa, nArray, h1Array, h2Array, boundary1IncludedArray, boundary2IncludedArray] = ...
    %               find_interval_overlaps2(f1Array, f2Array, boundary1IncludedArray, boundary2IncludedArray)
    %           PRO: Might not be so complicated. Can probably see every "interruption" b as b-, b, and b+, MAYBE.
    %   PROPOSAL: Always include boundaries in input intervals for overlap.
    %       NOTE: Can exclude boundaries from input intervals by excluding output intervals with length=0.
    
    % NOTE: Function is presently mathematically inconsistent (bad sign).
    %   Input:  Sets on form [a,b]
    %   Output: Sets on form (a,b) and [a,a], if output is interpreted as non-overlapping intervals.
    %   Ex: a<b<c
    %       [a,c], [b,b] --> [a,b), [b,b], (b,c]
    
    % Naming convention
    % =================
    % II = Input Interval
    
    PERMIT_ZLI_ARGS = 1;
    
    % Shorten variables.
    f1 = f1Array;
    f2 = f2Array;
    
    % ASSERTIONS
    assert(all(isfinite(f1)))
    assert(all(isfinite(f2)))
    assert(PERMIT_ZLI_ARGS || all(f1<f2))    % Optionally forbid input ZLIs.
    assert(all(f1<=f2), 'f1 <= f2 not true for all indices.')
    
    % Sorted list of "event times", i.e. when at least one interval begins or
    % ends.
    % NOTE: This eliminates any duplicate events, e.g. arising from identical
    %       input intervals.
    f12 = sort(unique([ f1(:); f2(:) ]));
    nF12 = numel(f12);
    
    % Sorted lists of beginning and end of output intervals.
    h1 = zeros(0,1);
    h2 = zeros(0,1);
    
    % "Current" set of (indices to) input intervals.
    set    = [];    % Current set
    % CA of sets of (indices to) input intervals.
    setsCa = cell(0,1);

    if nF12 >= 1
        f_prev = f12(1);    % Location of previous "event".
    end
    %===========================================================================
    % Iterate over "events" (unique beginnings/ends of input intervals), in
    % *increasing order*
    % --------------------------------------------------------------------------   
    
    % Algorithm works itself from low to high f values, updating a state as it
    % progresses over events.
    % When passing over an event, the state will either imply
    % (1) only adding              IIs to   the state,
    % (2) only removing            IIs from the state, or
    % (3) both adding and removing IIs from the state.
    % In a sense, (3) means that events are split into two: (a) Actions taken
    % just before the event, and (b) actions taction just after.
    
    % Add and remove from "set", and
    % store intermediate set value for every output interval.
    %===========================================================================
    for i = 1:nF12
        % =====================================================================
        % Find input intervals that (1) begin and (2) end, respectively, at the
        % current event.
        % ---------------------------------------------------------------------
        % NOTE: Does not alter the algorithm state yet.
        % =====================================================================
        % Intervals which begin at this event. / Set to add.
        set1 = [find(f12(i) == f1)]';
        % Intervals which end   at this event. / Set to remove.
        set2 = [find(f12(i) == f2)]';
        
        assert(~isempty(union(set1, set2)))
        
        %==================================
        % Store the state in return values
        %==================================
        if ~isempty(set1) && i ~= 1
            % CASE: There at least one II lower boundary at this event.
            %       It is not the first event.
            % ==> "set" will change (and it has been previously set).
            % ==> Store current set and interval for which it applies.
            
            % Act on "set" before the event.
            h1(end+1, 1)     = f_prev;
            h2(end+1, 1)     = f12(i);
            setsCa{end+1, 1} = set(:);
            f_prev           = f12(i);
        end

        %==========================================================
        % Update the algorithm state: Add input intervals to "set"
        %==========================================================
        assert(isempty(intersect(set, set1)))
        set = union(set, set1);    % NOTE: Sorts return value.
        % CASE: "set" corresponds to intervals covering or touching f12(i).

        %==================================
        % Store the state in return values
        %==================================
        if ~isempty(set2)
            % CASE: There is at least on II upper boundary at this event.
            % ==> "set" will change.
            % ==> Store current set and interval for which it applies.
            h1(end+1, 1)     = f_prev;
            h2(end+1, 1)     = f12(i);
            setsCa{end+1, 1} = set(:);
            f_prev           = f12(i);
        end

        %===============================================================
        % Update the algorithm state: Remove input intervals from "set"
        %===============================================================
        assert(all(ismember(set2, set)))
        set = setdiff(set, set2);    % NOTE: Sorts return value.

        % CASE: "set" corresponds to intervals covering after f12(i) (between
        %       "events").
    end
    
    nArray = cellfun(@numel, setsCa);
    
    assert(isempty(set))
    assert(iscolumn(h1))
    assert(iscolumn(h2))
    if nF12 > 0
        assert(all(h1 <= h2))
        %assert(sum(h2-h1) == f12(end)-f12(1))   % May fail due to rounding.
    end
    
    h1Array = h1;
    h2Array = h2;
end
