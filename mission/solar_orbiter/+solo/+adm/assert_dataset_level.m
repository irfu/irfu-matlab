%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-08-02
%
function assert_dataset_level(datasetLevel)
% PROPOSAL: Support L2S.
% PROPOSAL: Function for finding out whether level or not (i.e. without
%           exception).
% PROPOSAL: Constant for list of levels.

irf.assert.castring_regexp(datasetLevel, '(L1|L1R|L2|L3|HK|CAL)')
end
