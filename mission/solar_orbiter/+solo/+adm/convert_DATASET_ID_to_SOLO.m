%
% Convert a DATASET_ID from using ROC-SGSE to RODP/SOLO prefix.
% Ignores/permits DATASET_ID already using the SOLO_ prefix.
%
%
% Initially created 2020-03-03 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function datasetId = convert_DATASET_ID_to_SOLO(datasetId)
datasetId = regexprep(datasetId, '^ROC-SGSE_', 'SOLO_');
end
