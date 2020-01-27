%
% Split a string into consecutive parts, each one corresponding to a regexp match.
%
% The algorithm will try to match the beginning of str to the first regexp, then continue to match the remainder of the
% string for each successive regexp.
% NOTE: A regexp that does match an empty string (e.g. 'a*'), may return an empty substring. This is natural in this
% application, but is maybe not default regexp behaviour.
%
% The primary purpose of this function is to make it easy to parse a string syntax.
%
% NOTE: Not always usable. Ex: Matching with restrictive regexp (one or several) at the end and hacing regexp that
% permits ~arbitrary string in the beginning/middle. Arbitrary string will match until the end, preventing matching the
% ending regular expressions.
%
%
% ARGUMENTS
% =========
% str        : String
% regexpList : Cell array of strings, each one containing a regexp. ^ at the beginning of a regexp will be ignored.
%              NOTE: The sequence of regexes must match every single character in str.
%
%
% RETURN VALUE
% ============
% subStrList : Cell array of strings, each being a match for the corresponding string in regexpList.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-25
%
function subStrList = regexp_str_parts(str, regexpList)
% PROPOSAL: Better function name.
% 
% PROPOSAL: Somehow improve to prevent arbitrary strings in the middle.
%   PROPOSAL: Simultaneously begin searching from end.
%   PROPOSAL: Associate order of matching with each regexp. Must go from edges toward the middle (a regexp must be
%       applied to a string that is adjacent to at least one already matched substring.
%       CON: Still does not work for regexp surrounded by two ~arbitrary strings.
%           CON: Impossible problem to solve, even in principle(?).
%               PROPOSAL: Allow to search for matching substrings inside interior substring (substring not bounded by
%                         already made matches).

    subStrList = cell(numel(regexpList), 1);
    
    remStr = str;   % rem = remainder.
    for i = 1:numel(regexpList)
        % IMPLEMENTATION NOTE: Option "emptystring" is important. Can otherwise not distinguish between (1) no match, or
        % (2) matching empty string, e.g. for regexp ' *'. This important e.g. for matching unimportant whitespace in a
        % syntax.
        subStr = regexp(remStr, ['^', regexpList{i}], 'match', 'emptymatch');  
        if isempty(subStr)
            error('Could not match regexp "%s" to beginning of the remainder of the string, "%s".', regexpList{i}, remStr)
        end
        subStr = subStr{1};
        
        remStr = remStr(numel(subStr)+1:end);        
        subStrList{i} = subStr;
    end
    
    % ASSERTION
    if ~isempty(remStr)
        error('BICAS:regexp_str_parts:Assertion', 'Argument str="%s" does not match the submitted regular expressions.', str)
    end
end
