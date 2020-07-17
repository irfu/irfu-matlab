% Classify DATASET_ID.
%
% Returns easy-to-use flags to make it easy to implement different handling for different DATASET_IDs.
%
%
% NOTE: Only classifies BICAS' input & output datasets. Assertion for other.
% NOTE: Function deliberately ignores Skeleton/data version if part of the string.
% IMPLEMENTATION NOTE: Still recognizes old ROC-SGSE datasets since they may be found in global attribute
% DATASET_ID in old test files.
%
function C = classify_DATASET_ID(datasetId)
    % NEED: Be able to generalize nicely to datasets other than BICAS L1/L1R-->L2 processing: sweeps, currents, HK.
    % NEED: Be able to generalize nicely to BIAS L3 datasets.
    %   TODO-DEC: Counts as LFR? What if also uses TDS to derive L3 (should not happen)?
    %
    % PROPOSAL: Use regexp instead.
    %   PRO: Can more easily handle old ROC-SGSE datasets.
    %
    % PROPOSAL: Implement assertion on DATASET_ID via this function.
    %   Ex: bicas.assert_DATASET_ID
    %   Ex: bicas.swmode_defs.assert_DATASET_ID
    %   CON: Requires strict matching.
    %   PRO: Does not spread out the knowledge of DATASET_IDs.
    %   PROPOSAL: Flag for obsoleted DATASET_IDs that may be found in input datasets. Caller decides how to
    %       respond.
    % NEED?!: Some way of determining whether an obsoleted and current DATASET_ID are equivalent.
    %
    % PROPOSAL: Generalize to work for all DATASET_IDs (BICAS-related and not). Put outside BICAS.
    % PROPOSAL: Return whether SOLO or ROC-SGSE prefix.



    EJ_library.assert.castring(datasetId)
    
    % One flag per type of input/output voltage data.
    % IMPLEMENTATION NOTE: Avoiding the flag name isLfrCwf since it is ambiguous. isLfrSurvSwf is chosen in analogy with isLfrSurvCwf.
    C.isLfrSbm1    = 0;
    C.isLfrSbm2    = 0;
    C.isLfrSurvCwf = 0;
    C.isLfrSurvSwf = 0;
    C.isTdsCwf     = 0;
    C.isTdsRswf    = 0;
    % One flag per exact DATASET_ID.
    C.isCurrent    = 0;
    C.isBiasHk     = 0;
    % One flag per level.
    C.isL1         = 0;
    C.isL1R        = 0;
    C.isL2         = 0;



    % IMPLEMENTATION NOTE: Important to put L1R before L1. Otherwise, L1 will be matched to the beginning of L1R. ==>
    % Fail to match "R" (in L1R) to "_".
    [subStrList, remainingStr, perfectMatch] = EJ_library.str.regexp_str_parts(...
        datasetId, ...
        {'(SOLO|ROC-SGSE)', '_', '(HK|L1R|L1|L2)', '_', 'RPW-[A-Z12-]*'}, ...
        'permit non-match');
    
    % ASSERTION
    if ~perfectMatch
        % Better error message.
        error('Can not interpret datasetId="%s". The substring "%s" can not be interpreted', datasetId, remainingStr)
    end
    
    %prefix = subStrList{1};   % Currently not being used, but could potentially be.
    level  = subStrList{3};
    suffix = subStrList{5};
    
    switch(level)
        case 'L1'  ; C.isL1  = 1;
        case 'L1R' ; C.isL1R = 1;
        case 'L2'  ; C.isL2  = 1;
        case 'HK'  %; C.isHk  = 1;
        otherwise
            error('BICAS:proc_utils:Assertion:IllegalArgument', 'Can not handle DATASET_ID. datasetId="%s"', datasetId)
    end
    
    
    
    if     strcmp(datasetId, 'SOLO_L1_RPW-BIA-CURRENT')
        C.isCurrent = 1;
    elseif strcmp(datasetId, 'SOLO_HK_RPW-BIA')
        C.isBiasHk  = 1;
    else        
        if (C.isL1R | C.isL2)
            assert(strcmp(suffix(end-1:end), '-E'))
            suffixETruncated = suffix(1:end-2);
        else
            suffixETruncated = suffix;
        end
        
        switch(suffixETruncated)
            case 'RPW-LFR-SBM1-CWF' ; C.isLfrSbm1    = 1;
            case 'RPW-LFR-SBM2-CWF' ; C.isLfrSbm2    = 1;
            case 'RPW-LFR-SURV-CWF' ; C.isLfrSurvCwf = 1;
            case 'RPW-LFR-SURV-SWF' ; C.isLfrSurvSwf = 1;
            case 'RPW-TDS-LFM-CWF'  ; C.isTdsCwf     = 1;
            case 'RPW-TDS-LFM-RSWF' ; C.isTdsRswf    = 1;
            otherwise
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Can not handle DATASET_ID. datasetId="%s"', datasetId)
        end
    end

    %================================================
    % Set flags that can be derived from other flags
    %================================================
    C.isLfr = C.isLfrSbm1 | C.isLfrSbm2 | C.isLfrSurvCwf | C.isLfrSurvSwf;
    C.isTds = C.isTdsCwf  | C.isTdsRswf;
    C.isCwf = C.isLfrSbm1 | C.isLfrSbm2 | C.isLfrSurvCwf | C.isTdsCwf;
    C.isSwf = C.isLfrSurvSwf                             | C.isTdsRswf;

    
    
    % ASSERTION
    EJ_library.assert.struct(C, {...
        'isLfrSbm1', ...
        'isLfrSbm2', ...
        'isLfrSurvCwf', ...
        'isLfrSurvSwf', ...
        'isTdsCwf', ...
        'isTdsRswf', ...
        'isCurrent', ...
        'isBiasHk', ...
        'isL1', ...
        'isL1R', ...
        'isL2', ...
        'isCwf', ...
        'isSwf', ...
        'isLfr', ...
        'isTds'}, {})
end
