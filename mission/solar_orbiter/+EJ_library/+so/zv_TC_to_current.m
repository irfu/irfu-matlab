%
% Utility function for easily deriving (set) current at arbitrary timestamps from TC, i.e. almost directly from
% zVariables in SOLO_L1_RPW-BIA-CURRENT datasets. Therefore has some extra features to handle the quirks of the CDF
% format.
%
% Can be useful for plotting.
% NOTE: Used by BICAS, which want to detect duplicatesAnomaly and react to it itself.
%
% NOTE: Assumes that fill value is NaN. 
% NOTE: Ignores timestamps+bias when zvIBIASx1==NaN, because of the dataset format.
% NOTE: Tries to mitigate bug(?) in SOLO_L1_RPW-BIA-CURRENT that includes duplicate bias settings.
%
%
% ARGUMENTS
% ==========
% t1        : 1D vector. Time. Sorted. Type unimportant (linear increasing with time). Can be TT2000.
% zvIBIASx1 : 1D vector. Same length as t1. New bias values.
% t2        : 1D vector. Time. Same type of time as t1.
%
% 
% RETURN VALUE
% ============
% zvIBIASx2         : Current values at t2. Double. NaN for timestamps before first input timestamp.
% duplicatesAnomaly : Whether detected known anomaly in SOLO_L1_RPW-BIA-CURRENT datasets.
%                     RATIONALE: Caller (e.g. BICAS) can give error, warning.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
function [zvIBIASx2, duplicatesAnomaly] = zv_TC_to_current(t1, zvIBIASx1, t2)
    % PROPOSAL: Accept all three antennas at the same time.
    %   PRO: Can verify format better.
    % PROPOSAL: Accept exact CDF format.
    %   NOTE: Non-NaN fill value.
    %
    % PROPOSAL: Separate functions for converting
    %       ~zVar --> step function
    %       step function--> interpolation.
    
    
    
    % ASSERTIONS
    assert(numel(t1) == numel(zvIBIASx1), 'Arguments t1 and zvIBIASx1 do not have the same number of elements.')
    % NOTE: Do not check for MONOTONIC increase (yet) since it may be because of duplicate bias settings. Still useful
    % check since global.
    assert(issorted(t1), 'Argument t1 is not sorted.')
    
    % NOTE: (1) interp1 requires float or double, and (2) return value has to be float to store NaN anyway.
    t1        = double(t1);
    t2        = double(t2);
    zvIBIASx1 = double(zvIBIASx1);
    
    % Remove indices at which CURRENTS (not Epoch) are NOT NaN, i.e. which do not represent this antenna.
    % NOTE: Antenna is determined by the data in zvIBIASx1.
    bKeep     = ~isnan(zvIBIASx1);
    t1        = t1(bKeep);
    zvIBIASx1 = zvIBIASx1(bKeep);

    %======================================================================================================
    % CDF ASSERTION
    % Handle non-monotonically increasing Epoch
    % -----------------------------------------
    % NOTE: This handling is driven by
    % (1) wanting to check input data
    % (2) interp1 does not permit having identical x values/timestamps, not even with identical y values.
    %======================================================================================================
    if ~issorted(t1, 'strictascend')
        % CASE: Timestamps do NOT increase monotonically.
        
        % Set bIdent = wheter component is followed by identical value.
        bIdent = diff(t1) == 0;
        bIdent = [bIdent(:); false];   % Add last component to maintain same vector length.
        iIdent = find(bIdent);
        
        % ASSERTION: Identical timestamps correspond to identical bias settings.
        assert(all(zvIBIASx1(iIdent) == zvIBIASx1(iIdent+1)), ...
            'TC_to_current:Assertion', ...
            'Bias currents contain non-equal current values on equal timestamps on the same antenna.');
        
        % Mitigate: Remove duplicates
        t1        = t1(~bIdent);
        zvIBIASx1 = zvIBIASx1(~bIdent);
        duplicatesAnomaly = 1;
        
        % ASSERTION: Epoch increases monotonically (after mitigation)
        assert(issorted(t1, 'strictascend'), 'TC_to_current:Assertion', ...
            'Bias current timestamps do not increase monotonically after removing duplicate bias settings.')
    else
        % CASE: Timestamps do increase monotonically.
        duplicatesAnomaly = 0;
    end
    
    
    
    % IMPLEMENTATION NOTE: Bias currents are set VERY RARELY. Must therefore use interpolation method 'previous'.
    % NOTE: interp1 does NOT require t1 to be sorted.
    % NOTE: interp1 does not permit identical x1 values, not even if the coresponding y1 values are identical.
    zvIBIASx2 = interp1(t1, zvIBIASx1, t2, 'previous');
    
    % NOTE: interp1(... 'previous') returns NaN also for t2 > max(t1)!! Must therefore do this oneself.    
    zvIBIASx2(t2 > t1(end)) = zvIBIASx1(end);
end
