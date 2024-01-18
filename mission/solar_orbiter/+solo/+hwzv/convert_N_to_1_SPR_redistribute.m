%
% Convert zVariable-like variable from N samples/record to 1 sample/record (from
% a matrix to a column vector). Increases number of records.
%
% NOTE: This function is primarily intended to be used on zVariables.
%
%
% ARGUMENT
% ========
% zv1 : (iRecord, iSnapshotSample, iChannel)
%
%
% RETURN VALUE
% ============
% zv2 : (iRecord, iChannel).
%       Same number of components, nRecords2 = nRecords1*nSpr1
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-28, based on older code.
%
function zv2 = convert_N_to_1_SPR_redistribute(zv1)

irf.assert.sizes(zv1, [NaN, NaN, NaN])

zv  = permute(zv1, [2,1,3]);
zv2 = reshape(zv, size(zv,1) * size(zv,2), size(zv,3));

irf.assert.sizes(zv2, [NaN, NaN])
end
