%
% Only keep DSMDs with any of the specified DATASET_IDs.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-08.
%
function DsmdArray2 = filter_DSMD_DATASET_ID(DsmdArray1, datasetIdCa)
% PROPOSAL: Return indices instead.
%   CON: (Even more) trivial for the caller to do what function does.
% PROPOSAL: Convert to DSMD method.
%   PROPOSAL: Instance method?
%   PROPOSAL: Static method?

assert(isa(DsmdArray1, 'solo.adm.DSMD'))
irf.assert.castring_set(datasetIdCa)

bKeep = ismember({DsmdArray1.datasetId}, datasetIdCa);

DsmdArray2 = DsmdArray1(bKeep);
end
