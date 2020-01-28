function subsetStructList = select_cell_array_structs(structList, fieldName, selectionValuesList)
% Generic function. Given a cell array of structs, select the subset of these where a specified field matches any of the
% specified string values. 
%
% ASSUMES: Every value in the selection values list to be found exactly once (otherwise assertion error).
% ASSUMES: The structure field values are strings.
%
% NOTE: Can be used multiple times successively to select for two matching struct fields, e.g. datasetId
% AND skeletonVersionStr.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-19
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% structList          : Cell array of structs
% fieldName           : Name of a struct field.
% selectionValuesList : Cell array of strings.
% subsetStructList    : Cell array of structs. That subset of structs in "structList" where the field "fieldName"
%                       equals any of the string values in "selectionValuesList".



% Convert to cell array of strings.
structListValues = {};
for iStruct = 1:length(structList)
    structListValues{iStruct} = structList{iStruct}.(fieldName);
end

for iSelVal = 1:length(selectionValuesList)    
    jFound = find(strcmp(selectionValuesList{iSelVal}, structListValues));
    
    % ASSERTION
    if numel(jFound) ~= 1
        error('BICAS:select_cell_array_structs:Assertion:IllegalArgument', ...
            'Can not find exactly one structure with field "%s" with value "%s" (found %i of them).', ...
            fieldName, selectionValuesList{iSelVal}, numel(jFound))
    end
    
    subsetStructList{iSelVal} = structList{jFound};
end

end
