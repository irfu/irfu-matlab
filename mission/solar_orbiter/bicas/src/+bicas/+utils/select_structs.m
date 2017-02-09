% select_structs   Given a cell array of structs, select a subset by matching (string) fields.
% 
% Generic function. Given a cell array of structs, select the subset of these where a specified field matches any of the
% specified string values.
%
%
% NOTE: Can be used multiple times successively to select for two matching struct fields, e.g. dataset_ID
% AND skeleton_version_str.
%
% ASSUMES: All structure field values (in "structList") for the chosen field ("fieldName") are unique.
% ASSUMES: The structure field values are strings.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-19
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% structList       : Cell array of structs
% fieldName        : Name of a struct field.
% selectionValues  : Cell array of strings.
% subsetStructList : Cell array of structs. That subset of structs in "structList" whose field name "fieldName"
%                    equals any of the string values in "selectionValues".
%
function subsetStructList = select_structs(structList, fieldName, selectionValues)

% Convert to cell array of strings.
structListValues = {};
for iStruct = 1:length(structList)
    structListValues{iStruct} = structList{iStruct}.(fieldName);
end

for iSelVal = 1:length(selectionValues)    
    jFound = find(strcmp(selectionValues{iSelVal}, structListValues));
    
    % ASSERTION
    if numel(jFound) ~= 1
        error('select_structs:Assertion:IllegalArgument', 'Can not find exactly one structure with field "%s" with value "%s" (found %i of them).', ...
            fieldName, selectionValues{iSelVal}, numel(jFound))
    end
    
    subsetStructList{iSelVal} = structList{jFound};
end

end
