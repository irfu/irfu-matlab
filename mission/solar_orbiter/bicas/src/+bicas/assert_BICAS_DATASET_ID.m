%
% Assert DATASET_ID associated with BICAS, i.e. not just any DATASET_ID.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-08-02
%
function assert_BICAS_DATASET_ID(datasetId)
    % PROPOSAL: Implement via list.
    
    % '(ROC-SGSE|SOLO)_(L[12].+RPW-(LFR|TDS)-(SBM[12]|SURV|LFM)-(C|S|RS)WF.*)|HK_RPW_BIA)'
    % IMPLEMENTATION NOTE: Does not cover everything.
    % NOTE: MATLAB regxp can not handle recursive brackets it seems, i.e.
    % ((...|...)|(...|...)) .
    
    % NOTE: solo.adm.disassemble_DATASET_ID does some assertions on
    % sourceName and level. /2020-09-29
    [~, ~, descriptor] = solo.adm.disassemble_DATASET_ID(datasetId);
    
    % NOTE: Constrain DATASET_ID to roughly BICAS-related datasets.
    irf.assert.castring_regexp(descriptor, 'RPW-(BIA|LFR|TDS)[A-Z0-2-]*')
    
end
