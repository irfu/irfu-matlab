function s = merge_structs(s1, s2)
% Merge two structs by creating a new struct with the union of the fields and values.
% The field names must not overlap.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% s1, s2 : Two structs with non-overlapping lists of field names. Could be struct arrays (1D).
%          If one is a struct array, then the other one must be
%          either a scalar (1x1) struct, or have the same array size as the other struct.
% s      : Struct ROW array (regardless of sizes of s1, s2).
%
%
% NOTE: Can be used to effectively add a field to an existing struct array.
%
%
% Author: Erik P G Johansson, Sweden
% First created 2017-10-12

% "Test code":
% EJ_library.utils.merge_structs(struct('a', 1), struct('b', 2, 'c', 4))
% s1=struct('a', {1,2,3,4}); s2=struct('b', {2,3,4,6}', 'c', {4,5,6,7}'); s=EJ_library.utils.merge_structs(s1, s2)


    if ~isvector(s1) || ~isvector(s2)
        error('BICAS:merge_structs:Assertion', 'At least one of the structs is not a vector array.')
    end
    if ~(isscalar(s1) || isscalar(s2))
        if numel(s1) ~= numel(s2)
            error('BICAS:merge_structs:Assertion', 'Structs are arrays of different non-scalar size.')
        end
    end

    fieldNamesList1 = fieldnames(s1);
    structArgsList1 = {};
    for i = 1:length(fieldNamesList1)
        fieldName = fieldNamesList1{i};
        structArgsList1{end+1} = fieldName;
        structArgsList1{end+1} = {s1.(fieldName)};
    end
    
    fieldNamesList2 = fieldnames(s2);
    structArgsList2 = {};
    for i = 1:length(fieldNamesList2)
        fieldName = fieldNamesList2{i};
        structArgsList2{end+1} = fieldName;
        structArgsList2{end+1} = {s2.(fieldName)};
    end
    
    structArgsList = [structArgsList1, structArgsList2];    % Merge two row cell arrays.
    s = struct(structArgsList{:});    % Gives error if reuses same fieldname.
end
