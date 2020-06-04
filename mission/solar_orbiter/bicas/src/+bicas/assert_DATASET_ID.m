%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-08-02
%
function assert_DATASET_ID(datasetId)
% PROPOSAL: Optional second argument for pipeline ID.
% PROPOSAL: Implement via list.
% PROPOSAL: Implement via EJ_library.so.adm.classify_DATASET_ID. See implementation.

% '(ROC-SGSE|SOLO)_(L[12].+RPW-(LFR|TDS)-(SBM[12]|SURV|LFM)-(C|S|RS)WF.*)|HK_RPW_BIA)'
% IMPLEMENTATION NOTE: Does not cover everything. It is hard to cover "everything" MATLAB regxp can not
% handle recursive brackets it seems, i.e. ((...|...)|(...|...)) .
EJ_library.assert.castring_regexp(datasetId, '(ROC-SGSE|SOLO)_(L[12].?|HK)_RPW-(BIA|LFR|TDS)[A-Z1-2-]*')

end
