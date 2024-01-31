%
% Assert DSI associated with BICAS, i.e. not just any DSI.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-08-02
%
function assert_BICAS_DSI(dsi)
% PROPOSAL: Implement via list.

% '(ROC-SGSE|SOLO)_(L[12].+RPW-(LFR|TDS)-(SBM[12]|SURV|LFM)-(C|S|RS)WF.*)|HK_RPW_BIA)'
% IMPLEMENTATION NOTE: Does not cover everything.
% NOTE: MATLAB regxp can not handle recursive brackets it seems, i.e.
% ((...|...)|(...|...)) .

% NOTE: solo.adm.disassemble_DATASET_ID does some assertions on
% sourceName and level. /2020-09-29
[~, ~, descriptor] = solo.adm.disassemble_DATASET_ID(dsi);

% NOTE: Constrain DSI to roughly BICAS-related datasets.
irf.assert.castring_regexp(descriptor, 'RPW-(BIA|LFR|TDS)[A-Z0-2-]*')

end
