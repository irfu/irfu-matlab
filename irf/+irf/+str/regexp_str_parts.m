%
% Split a string into consecutive parts, each one corresponding to a regexp
% match.
%
% The primary purpose of this function is to make it easy to parse a string
% syntax.
%
%
% ALGORITHM
% =========
% The algorithm will try to match the beginning of the string to the first
% regexp, then continue to match the remainder of the string for each successive
% regexp. Perfect match means that the last regexp matches the then remaining
% string exactly.
% --
% NOTE: The matching can fail in two ways:
% (a) Algorithm runs out of string before running out of regular expressions.
% (b) Algorithm runs out of regular expressions before running out of string.
% --
% NOTE: A regexp that does match an empty string (e.g. 'a*'), may return an
% empty substring. This is natural in this application, but is maybe not default
% regexp behaviour.
% --
% NOTE: The algorithm does not work for "all" applications.
% Ex: Matching with restrictive regexp (one or several) at the end of string
% while simultaneously having regexp that permits ~arbitrary string in the
% beginning/middle. The arbitrary string will match until the end (maximal
% munch), preventing matching the last regular expressions.
%
%
% ARGUMENTS
% =========
% str
%       String
% regexpCa
%       Cell array of strings, each one containing a regexp. "^" at the
%       beginning of a regexp will be ignored.
%       NOTE: The sequence of regexes must match every single character in str.
% nonMatchPolicy
%       String constant determining what happens in the event of a non-perfect
%       match (including no match).
%           'assert match'
%           'permit non-match'
%       This refers to both kinds of failure (above).
%
%
% RETURN VALUE
% ============
% subStrCa
%       Cell array of strings, each being a match for the corresponding string
%       in regexpCa.
% remainingStr
%       The remainder of argument str that was not matched.
% isPerfectMatch
%       Logical. Whether matched all regular expressions to entire string.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-01-25.
%
function [subStrCa, remainingStr, isPerfectMatch] = regexp_str_parts(...
  str, regexpCa, nonMatchPolicy)

% PROPOSAL: Better function name.
%   PROPOSAL: ~token
%   PROPOSAL: read_tokens (plural)
%   PROPOSAL: read_regexp_tokens
%   PROPOSAL: read_arbitrary_tokens
%   PROPOSAL: tokenize
%       CON: Implies similarity between tokens, not arbitrary regular
%            expressions.
%
% PROPOSAL: Rename substrings "tokens".
%
% PROPOSAL: Somehow improve to prevent failing for regexps matching ~arbitrary strings before other regexp(s).
%   Ex: Matching SO dataset filename. Can not match "-CDAG" or ".cdf".
%   PROPOSAL: Option to search backwards (only).
%       TODO-DEC: How represent partial match (only some substrings)?
%   PROPOSAL: Simultaneously searching from both beginning and end.
%   PROPOSAL: Associate order of matching with each regexp. Must go from edges toward the middle (a regexp must be
%       applied to a string that is adjacent to at least one already matched substring.
%       CON: Still does not work for regexp surrounded by two ~arbitrary strings.
%           CON: Impossible problem to solve, even in principle(?).
%               PROPOSAL: Allow to search for matching substrings inside interior substring (substring not bounded by
%                         already made matches).
%                   NOTE: Not exact solution, but probably works practically for most applications.
%                   CON: Better done manually or by another function that uses this function. This is outside of
%                        this function's natural scope.
%   PROPOSAL: Two functions
%       [subStrCa, remainingStr] = regexp_str_parts(str, regexpCa, searchDirection, nonMatchPolicy)
%           Only search either forward or backward (argument).
%       [subStrCa, remainingStr] = regexp_str_parts2(str, regexpCa, searchPriorities, nonMatchPolicy)
%           Uses regexp_str_parts. Each regexp is associated with a number. Search forward or backward from the
%           edges depending on searchPriorities. searchPriorities = one integer per regexp. Can not be arbitrary;
%           must increase toward the middle.
%
% PROPOSAL: Do not throw error on not matching, but return special value.
%   PROPOSAL: Return empty.
%   PROPOSAL: Have policy argument determine whether exception or special return value.
%       PROPOSAL: 'assert match', 'permit no match'/'no assert match'
%           NOTE: Easy for caller to assert: assert(~isempty(subStrCa)).
%   PROPOSAL: Permit returning partial results and letting the caller know where it failed.
%
% PROPOSAL: Somehow specify whether to ignore case or not, for every regexp separately(!).
%
% PROPOSAL: Additional easy-to-use return value for verifying perfect match.
%   PROBLEM: How make it fit with assertionPolicy? Should not be in the way with 'assert match'.
%       PROPOSAL: [subStrCa, remainingStr, isPerfectMatch]
%           PRO: Backward-compatible.
%       PROPOSAL: [subStrCa, isPerfectMatch, remainingStr]
%           PRO: Does not need to store remainingStr, even if only wants isPerfectMatch.
%
% PROPOSAL: Implement using irf.str.read_token.
%
% NOTE: Could almost(?) use function to implement equivalent functionality of regular expressions with
% (positive) lookbehind+lookahead (other function).
%   PROPOSAL: Implement (other function) that uses regexp with effectively
%           (positive) lookbehind+lookahead to search through strings:
%       (1) Join three regexps (lookbehind, search pattern, lookahead) into one regexp, and search for matches.
%       (2) Use regexp_str_parts to find out which parts correspond to which of the three regexps.
%       (3) Return only the match for the search pattern regexp.
%       TODO-NI: Rigorous? Are there no cases it can not handle?!
%       TODO-NI: Can modify to use negative lookbehind+lookahead?



%==================================
% Interpret, verify nonMatchPolicy
%==================================
switch(nonMatchPolicy)
  case 'assert match'
    assertMatch = true;

  case 'permit non-match'
    assertMatch = false;

  otherwise
    % ASSERTION
    error('Illegal argument nonMatchPolicy="%s"', nonMatchPolicy)
end
clear nonMatchPolicy



%====================
% MATCHING ALGORITHM
%====================
subStrCa     = cell(0, 1);
remainingStr = str;
for i = 1:numel(regexpCa)
  % IMPLEMENTATION NOTE: Option "emptystring" is important. Can otherwise
  % not distinguish between (1) no match, or (2) matching empty string,
  % e.g. for regexp ' *'. This important e.g. for matching unimportant
  % whitespace in a syntax.
  subStrCell = regexp(remainingStr, ['^', regexpCa{i}], 'match', 'emptymatch');
  if isempty(subStrCell)
    % NOTE: isempty() REFERS TO CELL ARRAY returned from "regexp", not
    % any potential string. subStrCell should always be a cell array,
    % both for matches and non-matches, as well as for even empty string
    % match.

    % CASE: Failed to match regexp. NOTE: This includes the case of
    % "running out of string".

    if assertMatch
      % ASSERTION
      error('regexp_str_parts:Assertion', ...
        ['Could not match regular expression "%s" to', ...
        ' the beginning of the remainder of the string, "%s".'], ...
        regexpCa{i}, remainingStr)
    else
      % NOTE: subStrCa partially completed.
      isPerfectMatch = false;
      return
    end
  end
  % CASE: Successful match (has already left loop if not)
  subStr = subStrCell{1};

  remainingStr  = remainingStr(numel(subStr)+1:end);
  subStrCa{i,1} = subStr;    % Add to list of matches.
end
% CASE: Matched all regexps (but not necessarily against entire string).



%=======================================
% Check if algorithm matched everything
%=======================================
if ~isempty(remainingStr)
  % CASE: subStrCa "completed" (one string for every regexp), but
  % remainingStr is not empty.
  if assertMatch
    % ASSERTION
    error('regexp_str_parts:Assertion', ...
      ['Only the beginning of argument str="%s" matches', ...
      ' the submitted regular expressions.'], ...
      str)
  else
    isPerfectMatch = false;
    return
    % Do nothing
  end
end

isPerfectMatch = true;
end
