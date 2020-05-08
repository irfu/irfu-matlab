%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-08-02
%
function assert_dataset_level(datasetLevel)
% PROPOSAL: Have depend on ROC-SGSE/RODP pipeline?

EJ_library.assert.castring_regexp(datasetLevel, '(L2|L2S)')
end

