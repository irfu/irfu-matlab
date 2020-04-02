% Classify DATASET_ID.
%
% Returns easy-to-use flags to make it easy to implement different handling for different DATASET_IDs.
%
%
% NOTE: Does not recognize HK datasets.
% NOTE: Only classifies BICAS voltage input & output datasets. (Is there a good reason for this?)
% NOTE: Function deliberately ignores Skeleton_version.
% IMPLEMENTATION NOTE: Still recognizes old ROC-SGSE datasets since they may be found in global attribute
% DATASET_ID in old test files.
%
function C = classify_DATASET_ID(datasetId)
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
    
    
    
    EJ_library.assert.castring(datasetId)
    
    % One flag per type of input/output voltage data.
    % IMPLEMENTATION NOTE: Avoiding the flag name isLfrCwf since it is ambiguous. isLfrSurvSwf is chosen in analogy with isLfrSurvCwf.
    C.isLfrSbm1    = 0;
    C.isLfrSbm2    = 0;
    C.isLfrSurvCwf = 0;
    C.isLfrSurvSwf = 0;
    C.isTdsCwf     = 0;
    C.isTdsRswf    = 0;
    % One flag per level.
    C.isL1         = 0;
    C.isL1R        = 0;
    C.isL2         = 0;
    
    
    
    % IMPLEMENTATION NOTE: Important to put L1R before L1. Otherwise, L1 will be matched to L1R. ==> Fail to match R to
    % "_".
    [subStrList, remainingStr, perfectMatch] = EJ_library.utils.regexp_str_parts(...
        datasetId, ...
        {'(SOLO|ROC-SGSE)', '_', '(HK|L1R|L1|L2)', '_', 'RPW-[A-Z12-]*'}, ...
        'permit non-match');
    
    % ASSERTION
    if ~perfectMatch
        % Better error message.
        error('Can not interpret datasetId="%s". The substring "%s" can not be interpreted', datasetId, remainingStr)
    end
    
    %prefix = subStrList{1};   % Currently not being used, but could be.
    level  = subStrList{3};
    suffix = subStrList{5};
    
    switch(level)
        case 'L1'  ; C.isL1  = 1;
        case 'L1R' ; C.isL1R = 1;
        case 'L2'  ; C.isL2  = 1;
        otherwise
            error('BICAS:proc_utils:Assertion:IllegalArgument', 'Can not handle DATASET_ID. datasetId="%s"', datasetId)
    end
    
    if (C.isL1R | C.isL2)
        assert(strcmp(suffix(end-1:end), '-E'))
        suffix2 = suffix(1:end-2);
    else
        suffix2 = suffix;
    end
    
    switch(suffix2)
        case 'RPW-LFR-SBM1-CWF' ; C.isLfrSbm1    = 1;
        case 'RPW-LFR-SBM2-CWF' ; C.isLfrSbm2    = 1;
        case 'RPW-LFR-SURV-CWF' ; C.isLfrSurvCwf = 1;
        case 'RPW-LFR-SURV-SWF' ; C.isLfrSurvSwf = 1;
        case 'RPW-TDS-LFM-CWF'  ; C.isTdsCwf     = 1;
        case 'RPW-TDS-LFM-RSWF' ; C.isTdsRswf    = 1;
        otherwise
            error('BICAS:proc_utils:Assertion:IllegalArgument', 'Can not handle DATASET_ID. datasetId="%s"', datasetId)
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
        'isL1', ...
        'isL1R', ...
        'isL2', ...
        'isCwf', ...
        'isSwf', ...
        'isLfr', ...
        'isTds'}, {})
end
