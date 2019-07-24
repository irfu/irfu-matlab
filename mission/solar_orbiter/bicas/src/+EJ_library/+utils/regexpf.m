%
% Roughly like regexp, except
% (1) it only matches entire strings (not substrings). In practise, it surrounds the submitted regular expressions with ^ and $.
% (2) it returns a more usable logic array.
%
% regexpf = regexp (MATLAB builtin) + f (=full match)
% 
%
% ARGUMENTS
% ========= 
% regexPattern : CA string or cell array of CA strings. Each CA string is a regexp which may or may not be surrounded by
%                ^ and $.
% isMatch      : Logical array. True iff str is matched by corresponding regexp.
%                NOTE: This is different from "regexp" which returns a cell array of arrays.
%
%
% Initially created 2018-07-12 by Erik P G Johansson.
%
function isMatch = regexpf(str, regexPattern)
    
%    EJ_library.utils.assert.castring(str)
    
    if iscell(regexPattern)
        
        isMatch = cellfun(@(re) (~isempty(regexp(str, ['^', re, '$'], 'once'))), regexPattern);
        
    elseif ischar(regexPattern)
        
        %EJ_library.utils.assert.castring(regexPattern)
        
        isMatch = EJ_library.utils.regexpf(str, {regexPattern});    % RECURSIVE CALL
        
    end    
end
