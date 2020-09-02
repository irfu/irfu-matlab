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
    %
    %
    % PROPOSAL: Generalize to work for all DATASET_IDs (BICAS-related and not).
    % PROPOSAL: Return whether SOLO or ROC-SGSE prefix.
    % PROPOSAL: Function that splits up DATASET_ID string into parts:
    %   source name : SOLO, ROC-SGSE
    %   level       : L1 etc
    %   descriptor  : (HK) RPW-BIA etc
    %   TODO-DECISION: How should this relate to assertion on DATASET_ID?
    %
    %
    %
    % NOTE: In principle, this function is a substitute for multiple functions DATASET_ID-->boolean, which are defined
    % on only a subset of DATASET_IDs.
    % PROPOSAL: Have caller request flags. If a flag is not defined for the specified DATASET_ID, then assertion error.
    
    % PROPOSAL: Define constant table DATASET_ID-->set_of_values. Values can be
    %           false/true or other, or be undefined. Define static
    %           methods/functions that use table to translate table to relevant
    %           flag. If requested flag is undefined, then assertion error.
    %   PROPOSAL: Have in EJ_library.so.constants.
    %   PROPOSAL: Format: containers.Map: DATASET_ID --> struct with one flag
    %   (or variable) per field. Different DATASET_IDs can have different fields.
    %       CON: Duplicates for SOLO and ROC-SGSE.
    %           PROPOSAL: Level + ~descriptor --> struct
    %   PROPOSAL: Can not use patterns in the constants and DATASET_IDs.
    %
    % PROPOSAL: Have class with only static methods + constants table.
    %   PROPOSAL: Name DSI_classif, DATASET_ID_classif, DSI, DATASET_ID
    %   PROPOSAL: Move convert_DATASET_ID_to_SOLO to class.
    %       PROPOSAL: Rename to convert_ROCSGSE_to_SOLO.
    %   PROPOSAL: Not use an actual constant table, but a private static method that returns all available flags for a
    %             given DATASET_ID.
    %       PROPOSAL: Table DATASET_ID descriptor+level --> values (instead of from DATASET_ID).
    %       PROPOSAL: Table DATASET_ID descriptor       --> values (instead of from DATASET_ID, or level).
    %           NOTE: Returning level would be using separate functionality.
    %   PROPOSAL: Use list of DATASET_IDs for assertion on valid DATASET_ID.
    %
    % PROPOSAL: Split up into multiple functions which cover subsets of DATASET_IDs. Their results can then be combined.
    %   CON: Can not implement assertions for bad DATASET_IDs.
    %       PROPOSAL: Policy argument for assertion.
    %   PRO: Not all flags always make sense.
    %   CON: When a flag should be defined or not is not well-defined.
    %       Ex: isLfr, isTds may be well-defined for union(LFR,TDS), but should
    %            they be defined for HK too?
    %       Ex: Should isLfr == ~isTds?
    %       Ex: Should isCwf == ~isSwf?
    %   --
    %   NEED: bicas.proc, bicas.proc_sub:    Classify L1/L1R LFR/TDS datasets.
    %   NEED: bicas.swmode_defs:             Assert valid BICAS input/output DATASET_ID.
    %   NEED: bicas.get_master_CDF_filename: Assert valid BICAS output DATASET_ID
    %   NEED: EJ_library.so.psp2:       LFR CWF/SBM1/SBM2/SWF to set Rx.
    %                                        Potentially classify any BICAS-related dataset.
    %   NEED: parse_dataset_filename:        Distinguish SOLO & ROC-SGSE.



    EJ_library.assert.castring(datasetId)
    
    % One flag per type of input/output voltage data.
    % IMPLEMENTATION NOTE: Avoiding the flag name isLfrCwf since it is ambiguous. isLfrSurvSwf is chosen in analogy with isLfrSurvCwf.
    C.isLfrSbm1    = false;
    C.isLfrSbm2    = false;
    C.isLfrSurvCwf = false;
    C.isLfrSurvSwf = false;
    C.isTdsCwf     = false;
    C.isTdsRswf    = false;
    % One flag per exact DATASET_ID.
    %C.isCurrent    = false;   % Not used
    %C.isBiasHk     = false;   % Not used
    % One flag per level.
    C.isL1         = false;
    C.isL1R        = false;
    C.isL2         = false;
    C.isL3         = false;



    % IMPLEMENTATION NOTE: Important to put L1R before L1. Otherwise, L1 will be matched to the beginning of L1R. ==>
    % Fail to match "R" (in L1R) to "_".
    [subStrList, remainingStr, perfectMatch] = EJ_library.str.regexp_str_parts(...
        datasetId, ...
        {'(SOLO|ROC-SGSE)', '_', '(HK|L1R|L1|L2|L3)', '_', 'RPW-[A-Z12-]*'}, ...
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
        case 'L1'  ; C.isL1  = true;
        case 'L1R' ; C.isL1R = true;
        case 'L2'  ; C.isL2  = true;
        case 'L3'  ; C.isL3  = true;
        case 'HK'   % Do nothing. There is "isBiasHk" instead.
        otherwise
            error('BICAS:proc_utils:Assertion:IllegalArgument', 'Can not handle DATASET_ID. datasetId="%s"', datasetId)
    end
    
    
    
    if     strcmp(datasetId, 'SOLO_L1_RPW-BIA-CURRENT')
        %C.isCurrent = true;
    elseif strcmp(datasetId, 'SOLO_HK_RPW-BIA')
        %C.isBiasHk  = true;
    elseif any(strcmp(datasetId, {'SOLO_L3_RPW-BIA-EFIELD', 'SOLO_L3_RPW-BIA-SCPOT'}))
        % Do nothing.
    else        
        if (C.isL1R || C.isL2)
            assert(strcmp(suffix(end-1:end), '-E'))
            suffixETruncated = suffix(1:end-2);
        else
            suffixETruncated = suffix;
        end
        
        switch(suffixETruncated)
            case 'RPW-LFR-SBM1-CWF' ; C.isLfrSbm1    = true;
            case 'RPW-LFR-SBM2-CWF' ; C.isLfrSbm2    = true;
            case 'RPW-LFR-SURV-CWF' ; C.isLfrSurvCwf = true;
            case 'RPW-LFR-SURV-SWF' ; C.isLfrSurvSwf = true;
            case 'RPW-TDS-LFM-CWF'  ; C.isTdsCwf     = true;
            case 'RPW-TDS-LFM-RSWF' ; C.isTdsRswf    = true;
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
        'isL1', ...
        'isL1R', ...
        'isL2', ...
        'isL3', ...
        'isCwf', ...
        'isSwf', ...
        'isLfr', ...
        'isTds'}, {})
end
