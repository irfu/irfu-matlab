function [subsetStructArray, iMatching] = select_array_structs(structArray, varargin)
% Generic function. Given an array of structs, select the subset of these where ALL specified fields separately have to
% match ANY one of a list of alternatives ("selection values") for the specific fields.
%
% NOTE: This function is similar to select_cell_array_structs, but not really analogous since:
%   (1) It permits non-string selection values.
%   (2) It does not require there to be exactly one match for every selection value (can be zero, one, or several).
%   (3) It returns indices to the matching structs.
%   (4) It permits matching multiple fields (not just multiple selection values).
%  
%
% NOTE: Not meant to be fast since it does not use fancy indexing code.
% NOTE: Uses "isequaln" for comparing field values and selection values. ==> NaN equals itself, 0==false, 1==true.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-03-15
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% structArray         : Array of structs.
% varargin            : A list of pairs of variables:
%     fieldName           : Name of a struct field.
%     selectionValuesList : Cell array of values to search for.
% subsetStructArray   : Array of structs. That subset of structs in "structArray" where the field "fieldName"
%                       equals any of the values in "selectionValuesList".
% iMatching           : Logical column vector the same length as structArray has elements. True for every index which
%                       matches. This is useful for modifying the original array of structs.
%

% ASSERTION
if mod(numel(varargin), 2) ~= 0
    error('BICAS:select_array_structs:Assertion', 'varargin does not consist of pairs of arguments.')
end



selectionArgumentList = varargin;
iMatchingSoFar = true(numel(structArray), 1);
while numel(selectionArgumentList) > 0
    
    fieldName           = selectionArgumentList{1};
    selectionValuesList = selectionArgumentList{2};
    if ~iscell(selectionValuesList)
        error('BICAS:select_array_structs:Assertion:IllegalArgument', 'selectionValuesList is not a cell array.')
    end
    
    iMatchingField = false(numel(structArray), 1);
    for iSelVal = 1:numel(selectionValuesList)
        
        % NOTE: Only iterate over those indices which match (so far) -- Faster.
        % IMPLEMENTATION NOTE: Must iterate over row vector of values. Otherwise the for statement will behave badly
        %                      for the case of zero-length vector.
        for iStruct = find(iMatchingSoFar)'
            if isequaln(...
                    structArray( iStruct ).( fieldName ), ...
                    selectionValuesList{ iSelVal }...
                    )
                iMatchingField(iStruct) = true;
            end
        end
    end
    
    % Effectively remove (set to false) those indices which do not match the last test field.
    iMatchingSoFar = and(iMatchingSoFar, iMatchingField);   

    selectionArgumentList = selectionArgumentList(3:end);
end



% Assign return values.
subsetStructArray = structArray(iMatchingSoFar);
iMatching         = iMatchingSoFar;
end





% 
% function [subsetStructArray, i] = select_array_structs(structArray, fieldName, selectionValuesList)
% % Generic function. Given an array of structs, select the subset of these where a specified field matches any of the
% % specified string values. 
% %
% % NOTE: This function is NOT entirely analogous to select_cell_array_structs since:
% %   (1) It permits non-string selection values
% %   (2) It does not require there to be exactly one match for every selection value (can be zero, one, or several).
% %
% % NOTE: Not meant to be fast since it does not use fancy indexing code.
% % NOTE: Uses "isequaln" for comparing field values and selection values. ==> NaN equals itself.
% %
% %
% % Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% % First created 2017-03-15
% %
% %
% % ARGUMENTS AND RETURN VALUES
% % ===========================
% % structArray         : Array of structs.
% % fieldName           : Name of a struct field.
% % selectionValuesList : Cell array of values to search for.
% % subsetStructArray   : Array of structs. That subset of structs in "structArray" where the field "fieldName"
% %                       equals any of the values in "selectionValuesList".
% 
% 
% 
% % Convert relevant field values to cell array.
% fieldValues = {};
% for iStruct = 1:length(structArray)
%     fieldValues{iStruct} = structArray(iStruct).(fieldName);
% end
% 
% iSave = [];
% for iSelVal = 1:length(selectionValuesList)
%     for iStruct = 1:length(structArray)
%         if isequaln(fieldValues{iStruct}, selectionValuesList{iSelVal})
%             iSave(end+1) = iStruct;
%         end
%     end
% end
% 
% subsetStructArray = structArray(iSave);
% end
