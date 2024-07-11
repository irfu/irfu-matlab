%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-08-02
%
function assert_dataset_level(datasetLevel)
irf.assert.castring_regexp(datasetLevel, '(L1|L1R|L2|L3|HK|CAL)')
end
