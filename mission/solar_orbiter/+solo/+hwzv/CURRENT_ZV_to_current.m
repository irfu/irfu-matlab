%
% Utility function for converting
%   from  ___ALMOST___ SOLO_L1_RPW-BIA-CURRENT dataset zVariables
%   to    easier-to-use arrays.
%
% Also checks the format, and checks for a known but "mitigatable" "duplicate
% bias anomaly" (see Section below). Therefore has some extra features to handle
% the quirks of the SOLO_L1_RPW-BIA-CURRENT dataset format.
%
% NOTE: "CURRENT" in the function name refers to SOLO_L1_RPW-BIA-CURRENT datasets.
% NOTE: Does not try to calibrate or convert units. The output bias is of the
% same type as the bias input.
% NOTE: Does not interpolate values to new timestamps.
% NOTE: See CURRENT_ZV_to_current_interpolate.
%
%
% DUPLICATE BIAS ANOMALY
% ======================
% There has historically been a mitigatable data anomaly in the CURRENT datasets
% in the form of repeated bias current values (same timestamp, same bias
% current). This was due to a ROC GitLab BICAS issue #17,
% https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/17
% The ROC issue should be fixed now. /Erik P G Johansson, 2020-09-15
%
%
% ALGORITHM / EFFECT
% ==================
% (1) Removes timestamp+bias when zvIBIASx1==NaN, since these legitimately occur
%     in the dataset format (it means that bias current is set on another
%     antenna in the same CDF record).
% (2) Removes duplicate bias settings (successive non-NaN data points with the
%     same timestamp and bias setting). Returns a flag value for whether this
%     has happened or not. This kind of data can be found due to the "duplicate
%     bias anomaly" (see Section) SOLO_L1_RPW-BIA-CURRENT and is therefore
%     flagged rather than asserted not to happen so that the caller can select
%     whether to give error/warn/accept, whether to use this as mitigation or
%     not.
% --
% NOTE: Does NOT remove bias settings that set the bias to the preceding value
% at a later timestamp ("unnecessary later bias settings").
%
%
% ARGUMENTS
% =========
% t1
%       Nx1 vector. Increasing (sorted; assertion), not necessarily
%       strictly. Time. Time representation unimportant as long as
%       increases with time. Can be e.g. TT2000.
% zvIBIASx1
%       Nx1 vector. Floating-point. Same length as t1. Bias values.
%       NOTE: NaN is fill value (not e.g. -1e31).
%
%
% RETURN VALUES
% =============
% t2
%       1D vector. Time. Same type of time as t1.
%       Subset of timestamps t1. See algorithm.
% zvIBIASx2
%       1D vector. Current values at t2. Double.
% duplicatesAnomaly
%       Whether has detected known anomaly mentioned above in
%       SOLO_L1_RPW-BIA-CURRENT datasets. Iff 1/true, then found duplicate
%       timestamps, with the SAME bias current.
%       NOTE: Bias current is still unambiguous in this case and the
%       function is designed to handle this case.
%       RATIONALE: Caller (e.g. BICAS) can give error, warning.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-23.
%
function [t2, zvIBIASx2, duplicatesAnomaly] = CURRENT_ZV_to_current(t1, zvIBIASx1)

    % ASSERTIONS
    irf.assert.sizes(...
        t1,        [-1, 1], ...
        zvIBIASx1, [-1, 1]);
    % Require current datatype that can store NaN.
    assert(isfloat(zvIBIASx1))
    % IMPLEMENTATION NOTE: Do not check for STRICT increase (yet) since it might
    % not be so because of duplicate bias settings anomaly. Checking for
    % NON-STRICT incrementation is still a useful check since it is global.
    assert(issorted(t1), 'Argument t1 does not increase, is not sorted.')

    % NOTE: return value has to be float to store NaN anyway.
%     t1        = double(t1);
%     zvIBIASx1 = double(zvIBIASx1);

    % Remove indices at which CURRENTS (not Epoch) are NOT NaN, i.e. which
    % provide actual bias values on this antenna.
    % NOTE: Antenna is determined by the data in zvIBIASx1.
    % NOTE: Need to specify row index to ensure correct size for empty
    % result variables.
    bKeep     = ~isnan(zvIBIASx1);
    t1        = t1       (bKeep, 1);
    zvIBIASx1 = zvIBIASx1(bKeep, 1);

    %============================================================================
    % CDF ASSERTION
    % Handle non-strictly increasing Epoch
    % ------------------------------------
    % NOTE: This handling is driven by
    % (1) wanting to check input data
    % (2) interp1 does not permit having identical x values/timestamps, not even
    %     with identical y values.
    %============================================================================
    if ~issorted(t1, 'strictascend')
        % CASE: Timestamps do NOT increase strictly.

        % Set bDupl = whether component (timestamp) is followed by identical
        % value (duplicate).
        bDupl = (diff(t1) == 0);
        % Add last component to maintain same vector length.
        bDupl = [bDupl(:); false];
        iDupl = find(bDupl);

        % ASSERTION: Successive duplicate timestamps correspond to identical
        % bias settings.
        assert(all(zvIBIASx1(iDupl) == zvIBIASx1(iDupl+1)), ...
            'TC_to_current:Assertion', ...
            ['Bias currents contain non-equal current values on equal', ...
            ' timestamps on the same antenna.']);

        %=============================
        % Mitigate: Remove duplicates
        %=============================
        t1        = t1(~bDupl);
        zvIBIASx1 = zvIBIASx1(~bDupl);
        duplicatesAnomaly = 1;

        % ASSERTION: Epoch strictly increases (after mitigation)
        assert(issorted(t1, 'strictascend'), ...
            'CURRENT_ZV_to_current:Assertion', ...
            ['Bias current timestamps do not strictly increase after', ...
            ' removing duplicate bias settings.'])
    else
        % CASE: Timestamps do increase strictly.
        duplicatesAnomaly = 0;
    end

    t2        = t1;
    zvIBIASx2 = zvIBIASx1;
end
