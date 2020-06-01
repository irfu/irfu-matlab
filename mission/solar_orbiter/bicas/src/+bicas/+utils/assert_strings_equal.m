%
% Assert that the strings in stringList (cell array of strings) are equal and that there is at least one string.
% Zero strings count as assertion false (error/warning).
%
% NOTE: Can compare an arbitrary number of strings, not just two.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created <=2016-10-28
%
function assert_strings_equal( L, giveError, stringList, msg )
    % PROPOSAL: Abolish?
    
    uniqueStrings = unique(stringList);
    
    if numel(uniqueStrings) ~= 1   % NOTE: stringList empty ==> uniqueStrings empty ==> Warning/error
        
        if numel(stringList) == 0
            stringListString = '(no strings)';
        else
            stringListString = ['"', strjoin(stringList, '", "'), '"'];
        end
        
        bicas.utils.react_to_false_assertion(L, giveError, [msg, ' String values: ', stringListString])
    end
    
end
