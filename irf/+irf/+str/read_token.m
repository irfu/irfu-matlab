%
% Utility function for parsing char string, one regular expression-defined token
% at a time from beginning or end.
%
% Useful for complex parsing of strings.
% ** Ex: Parsing filenames
% ** Ex: Complete parsing of more complex text file syntaxes, e.g. ODL, CSV
%        tables with both quoted and unquoted values.
%
%
% ALGORITHM
% =========
% Find the first regexp in (ordered) list of regexps that matches the
% beginning/end of the string. Return substring succeeding/preceding the match.
% --
% DESIGN INTENT: The code deliberately does NOT throw exceptions when failing to
% match regexp so that the caller can throw customized error messages with
% location etc, or taking other non-error actions.
%
%
% ARGUMENTS
% =========
% str       : String
% searchDir :  1 : Search forward
%           : -1 : Search backward
% varargin  : Ordered list of regular expressions.
%             NOTE: An algorithm for parsing of text often requires a special
%             case for when argument "str" is empty. Note that this is
%             distinctly different from a regexp that COULD match an empty
%             substring, but DOES NOT HAVE TO cover the entire string, e.g.
%             regexp '.*' or 'X*'  (if it was not for maximal munch), or
%             '(|word)'. This can be implemented by using the regexp "$", at the
%             appropriate location in the list of regexps, and then handle that
%             in analogy with other regexp, e.g. by using switch(iRegexp).
%
%
% RETURN VALUES
% =============
% token        : Substring from the beginning of str (argument) that matches
%                regexp varargin{iRegexp}.
% remainingStr : The remainder of argument str.
% iRegexp      : Index into varargin for the matching regex.
%                -1 if no match was found.
% --
% NOTE: searchDir== 1 ==> [token, remainingStr] == str
% NOTE: searchDir==-1 ==> [remainingStr, token] == str
%
%
% Initially created 2020-03-10, by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function [token, remainingStr, iRegexp] = read_token(str, searchDir, varargin)
% PROPOSAL: Use for irf.PDS_utils.convert_ODL_to_structs.
%   NOTE: Function is a modified copy of nested function.
%
% PROPOSAL: Do not return remaining string, but index to where it begins.
%   PRO: Potentially faster.
%       CON: Various experience shows that indexing can be slow.
%   CON: Caller has to use index string to get remaining string, which it probably wants anyway.
%       PRO: Easier to give good error messages when has exact location of error in original string.
%           CON: Probably used to parse multiple tokens which it makes it more complicated to keep track of the
%                exact index in the original string. This would not help much.
%
% PROPOSAL: Be strict about empty strings always size 1x0.
%   PROPOSAL: General function for normalizing strings (forze empty strings to be 1x0).
%       PROPOSAL: normalize_str_1x0, empty_str_to_1x0, norm_str_1x0.
%
% PROPOSAL: Return values leftStr+rightStr (not token+remainingStr), regardless of search direction.
%   PRO: Easier to intuitively interpret return values (merge the into one string in one's head).
%   PRO: Easier to write assertion leftStr+rightStr==str
%   CON: Inconsistent for caller.
%       CON: Inconsistent in algorithm (does not come naturally).
%
% PROPOSAL: String constants (instead of numbers) for searchDir.
%



ES = char(zeros(1,0));

% Configuration based on search direction
% ---------------------------------------
% (1) Prefix/suffix regex with ^ or $.
% (2) How to extract remaining string.
switch(searchDir)
  case  1
    regexpModifFunc = @(regexpStr)  (['^', regexpStr]);
    remainingStr    = @(str, token) (str(numel(token)+1 : end));
  case -1
    regexpModifFunc = @(regexpStr)  ([regexpStr, '$']);
    remainingStr    = @(str, token) (str(1 : end-numel(token)));
  otherwise
    % ASSERTION
    error('Illegal argument searchDir="%g".', searchDir)
end



for iRegexp = 1:length(varargin)

  % IMPLEMENTATION NOTE: Want to permit regular expressions that match an
  % EMPTY string to do so. Therefore HAS TO USE the regexp-returned
  % starting/ending index of the match. The regexp-returned STRING is
  % empty for both
  % (1) matching empty string, and
  % (2) no match.
  %   NOTE: Unimportant whether uses
  % IMPLEMENTATION NOTE: Need 'emptymatch' for SOME empty matches to be
  % recognized
  [token, iMatch] = regexp(str, regexpModifFunc(varargin{iRegexp}), ...
    'match', 'start', 'once', 'emptymatch');

  if ~isempty(iMatch)
    % CASE: Found match.

    % IMPLEMENTATION NOTE: Needed for returned empty token to always be
    % size 1x0. Will otherwise only sometimes be.
    if isempty(token)
      token = ES;
    end

    remainingStr = remainingStr(str, token);

    % Can fail because of size 0x0 vs 1x0 empty strings.
    %assert(strcmp([token, remainingStr], str))

    return    % NOTE: Exit function.
  end
end
% CASE: Found no match.

token        = ES;
remainingStr = str;   % NOTE: If empty, then might not be size 1x0.
iRegexp      = -1;

end
