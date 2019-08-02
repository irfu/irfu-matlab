%
% Construct a DVID (DATASET_ID+skeleton version, ID) derived from a dataset ID and skeleton version
%   V01 + ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E
% --> V01_ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-08-02
%
function dvid = construct_DVID(datasetId, skeletonVersionStr)
% PROPOSAL: Move to ~constants (collect decision functions).

bicas.assert_DATASET_ID(datasetId)
bicas.assert_skeleton_version(skeletonVersionStr)

% IMPLEMENTATION NOTE: Has to work sensibly for both ROC-SGSE and RODP/SOLO dataset IDs.
% IMPLEMENTATION NOTE: Put skeleton version at beginning of PDID since PDIDs then line up better when
% printed in a list above each other. Easier to read.
dvid = sprintf('V%s_%s', skeletonVersionStr, datasetId);

end
