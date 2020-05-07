%
% Roughly like MATLAB's "regexp", except
% (1) it only matches entire strings (not substrings). In practice, it surrounds the submitted regular expressions with ^ and $.
% (2) it can match empty strings (option in regexp),
% (3) it returns a more usable logic array.
% (4) it can iterate over both strings, and regular expressions (but not simultaneously).
%
% regexpf = regexp (MATLAB builtin) + f (=full match)
% 
%
% ARGUMENTS
% ========= 
% str          : String or cell array of strings. Empty string must be of class char, i.e. e.g. not [] (due to function "regexp").
%                NOTE: Will permit empty strings to match a regular expression.
% regexPattern : String or cell array of strings. Each string is a regexp which may or may not be surrounded by
%                ^ and $.
% NOTE: str and regexPattern must not both simultaneously be non-scalar cell arrays.
%
%
% RETURN VALUE
% ============
% isMatch      : Logical array. True iff str is matched by corresponding regexp.
%                Same size and indices as any non-scalar argument cell array.
%                NOTE: This is different from "regexp" which returns a cell array of arrays.
%
%
% Initially created 2018-07-12 by Erik P G Johansson.
%
function isMatch = regexpf(str, regexpPattern)
% PROPOSAL: Allow both str and regexpPattern to be cell arrays of strings.
%   CON: Using both simultaneously means having to return a 2D matrix and the caller having to keep track of indices with different meanings.
%       PRO: Bad for backwards compatibility.
%           CON: So far, only regexpPattern has been allowed to be a cell array. ==> isMatch is always a 1D array,
%                and the callers probably do not care about the index (presumaby; verify first).
%       CON: Breaks existing functionality. regexpPattern is allowed to be an arbitrary dimension cell array (not just
%            1D array) and the return result should have the same dimension.
%           CON-PROPOSAL: Allow both str and regexpPattern to be arbitrary dimension cell arrays, as long as dimensions
%               do not conflict (must not have any overlapping non-size-one dimensions/indices).
%               CON/PROBLEM: Difficult to implement.


    
    % IMPLEMENTATION NOTE: Empty string must be of class char, i.e. not [] (due to function "regexp").
    % IMPLEMENTATION NOTE: regexp accepts cell array of strings for strings to match (not regexp) which this function is
    % not supposed to handle (at least not yet). Must therefore explicitly check that str is not a cell array.
    %EJ_library.assert.castring(str)
    
    %====================================================
    % Normalize input: Turn strings into 1x1 cell arrays
    %====================================================
    if ischar(str)
        str = {str};
    end
    if ischar(regexpPattern)
        regexpPattern = {regexpPattern};
    end
    
    
    
    if     numel(str) >= 1 && numel(regexpPattern) == 1
        % NOTE: Handles case of both str and regexpPattern being scalar.
        
        isMatch = cellfun(@(s) (~isempty(regexp(s, ['^', regexpPattern{1}, '$'], 'once', 'emptymatch'))), str);
        
    elseif numel(str) == 1 && numel(regexpPattern) > 1
    
        % IMPLEMENTATION NOTE: regexp option "emptymatch" is required to match empty strings.
        isMatch = cellfun(@(re) (~isempty(regexp(str{1}, ['^', re, '$'], 'once', 'emptymatch'))), regexpPattern);
    else
        % ASSERTION
        error('Both arguments are non-scalar (more than one element) cell arrays. Function can not handle that case.')
    end
    
    assert(islogical(isMatch))
end
