%
% Merge two structs by creating a new struct with the union of the fields and
% values. The field names must not overlap. Not recursive.
%
% NOTE: Implemented by creating one combined call struct(fieldName1,
% fieldValue1, ...).
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% S1, S2
%       Two structs with non-overlapping lists of field names. Could be struct
%       arrays (1D). If one is a struct array, then the other one must be either
%       a scalar (1x1) struct, or have the same array size as the other struct.
% S12
%       Struct ROW array (regardless of sizes of s1, s2).
%
%
% NOTE: Can be used to effectively add a field to an existing struct array.
% NOTE: irf.ds.add_struct_to_struct is similar but not the same.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-10-12
%
function S12 = merge_structs(S1, S2)
%
% PROPOSAL: Replace (reimplement?) with irf.ds.add_struct_to_struct.
% PROPOSAL: ATEST code.
% PROPOSAL: Make work with any array size.
%
% "Test code":
% bicas.utils.merge_structs(struct('a', 1), struct('b', 2, 'c', 4))
% s1=struct('a', {1,2,3,4}); s2=struct('b', {2,3,4,6}', 'c', {4,5,6,7}'); s=bicas.utils.merge_structs(s1, s2)


% ASSERTIONS
if ~isvector(S1) || ~isvector(S2)
  error('merge_structs:Assertion', ...
    'At least one of the structs is not a vector array.')
end
if ~(isscalar(S1) || isscalar(S2))
  if numel(S1) ~= numel(S2)    % NOTE: Only works since S1,S2 are vectors.
    error('merge_structs:Assertion', ...
      'Structs are arrays of different non-scalar size.')
  end
end

fieldNamesList1 = fieldnames(S1);
structArgsList1 = {};
for i = 1:length(fieldNamesList1)
  fieldName = fieldNamesList1{i};
  structArgsList1{end+1} = fieldName;
  structArgsList1{end+1} = {S1.(fieldName)};
end

fieldNamesList2 = fieldnames(S2);
structArgsList2 = {};
for i = 1:length(fieldNamesList2)
  fieldName = fieldNamesList2{i};
  structArgsList2{end+1} = fieldName;
  structArgsList2{end+1} = {S2.(fieldName)};
end

% Merge two row cell arrays.
structArgsList = [structArgsList1, structArgsList2];
% NOTE: Gives error if reuses same fieldname, i.e. fieldname collisions in
% S1 and S2.
S12 = struct(structArgsList{:});
end
