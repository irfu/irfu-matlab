%
% Function for "translating" a string into another value from a table. Give proper error message if no match.
%
%
% RATIONALE
% =========
% Primarily intended as a utility function to avoid common verbose switch-case statements (with an "otherwise"
% assertion) which interprets and verifies string constants, only to assign new values to some other variable(s) in
% every case statement. Can then write the code more in the form of a table. Can be used to assign multiple variables in
% every case by having e.g. cell array after values.
% 
%
% ARGUMENTS
% =========
% beforeAfterTable : Cell array of (a) cell arrays of strings, and () arbitrary values.
%                    {iMeaning, 1} = Cell array of strings. Is allowed to be empty but will then never match.
%                    {iMeaning, 2} = Arbitrary value to be returned.
%                    NOTE: One (probably) does not want there to be any duplicate strings.
%                    RATIONALE: Argument has this structure to make before+after clear when hardcoding it.
% beforeStr        : String.
% afterValue       : beforeAfterTable{iMeaning, 2} for which beforeAfterTable{iMeaning, 1}==beforeStr.
% NOTE: Empty string matches empty string.
%
%
% Initially created 2019-09-18 by Erik P G Johansson.
%
function afterValue = translate(beforeAfterTable, beforeStr, errorMsgId, errorMsg)
% PROPOSAL: Submit function returning error message string. Only evaluated if error.
% PROPOSAL: After value for no match.
% PROPOSAL: errorMsgId == '' ==> errorMsg is afterValue for match.
% PROPOSAL: One argument fewer ==> Last argument is afterValue if no match.
% PROPOSAL: Be able to submit separate before and after tables. ==> Varying number of arguments.
% PROPOSAL: Use functions as after values. Only evaluate if returned.
%   CON: Can effectively be used so already if the caller immediately evaluates the after value.
% PROPOSITION: Function is unnecessary.
%   PRO: Does not shorten the code enough to warrant a function.
%       PRO: Can write case statements on one row.
% PROPOSAL: Be able to use non-string beforeStr.
%   PRO: Can not always do with switch-case. ==> Can avoid if-elseif-elseif-...-else statements
%       Ex: Numeric vectors

[beforeTable, afterTable] = convert_before_after_table(beforeAfterTable);

assert(~isempty(errorMsgId))   % NOTE: Empty errorMsgId ==> error() will not throw exception.



for iMeaning = 1:numel(beforeTable)
    
    beforeAlts = beforeTable{iMeaning};
    
    % IMPLEMENTATION NOTE: ismember does not work as expected if beforeAlts is not a cell array, in particular a plain
    % string.
    assert(iscell(beforeAlts))
    
    if ismember(beforeStr, beforeAlts)
        afterValue = afterTable{iMeaning};
        return
    end
end

error(errorMsgId, errorMsg)

end



function [beforeTable, afterTable] = convert_before_after_table(beforeAfterTable)
assert(size(beforeAfterTable, 2) == 2)

beforeTable = beforeAfterTable(:,1);
afterTable  = beforeAfterTable(:,2);
end
