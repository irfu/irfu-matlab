%
% Function for "translating" a string into another value (not necessarily
% string) from a table. Give proper error message if no match.
%
%
% RATIONALE
% =========
% Primarily intended as a utility function to avoid common verbose switch-case
% statements (with an "otherwise" assertion) which interprets and verifies
% string constants, only to assign new values to some other variable(s) in every
% case statement. One can then write the code more in the form of a table.
% ** Can be used to assign multiple variables (output) in every case by having
%    e.g. cell array as values.
% ** Can use the same table structure to iterate over all "switch-case"
%    statements.
%
%
% ARGUMENTS
% =========
% SYNTAX: translate(... , errorMsgId, errorMsg)
%   No match will lead to error. (There can not be multiple matches.)
% SYNTAX: translate(... , nonMatchValue)
%   No match will be accepted.
% --
% table
%       Cell array of (a) cell arrays of strings, and (b) arbitrary values.
%       {iRule, 1} = Cell array of unique strings, "keys". Is allowed to be
%                    empty but will then never match.
%       {iRule, 2} = Arbitrary value to be returned.
%       NOTE: Asserts there are only unique keys across all rules.
%       RATIONALE: Argument has this structure to make keys+values clear when
%       hardcoding it using literals.
% key
%       String
% varargin
%       Alt 1: Permit non-match.
%         varargin{1} = nonMatchValue
%       Alt 2: Error on non-match.
%         varargin{1} = errorMsgId
%         varargin{2} = errorMsg
% --
% NOTE: Empty string matches empty string.
% NOTE: Counts '' and char(zeros(1,0)) as identical, both for matching and
% asserting unique keys (due to the behaviour of "ismember" and "unique", as
% opposed to "strcmp").
%
%
% RETURN VALUE
% ============
% value
%       Same as table{i, 2} for which table{i, 1}==key (string comparison).
%
%
% Initially created 2019-09-18 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function value = translate(table, key, varargin)
% PROPOSAL: Submit function returning error message string. Only evaluated if error.
%   PRO: Useful for complex error messages.
%
% PROPOSAL: Use functions as values. Only evaluate if returned.
%   CON: Can effectively be used so already if the caller immediately evaluates the value.
%
% PROPOSAL: Value for no match.
%   PROPOSAL: One argument fewer ==> Last argument is value if no match.
%
% PROPOSITION: Function is unnecessary.
%   PRO: Does not shorten the code enough to warrant a function.
%       PRO: Can write case statements on one row.
%   PRO: Function is only used by bicas.sprint_BSO().  /2024-07-23
%
% PROPOSAL: Be able to use non-string key.
%   PRO: Can not always do with switch-case. ==> Can avoid if-elseif-elseif-...-else statements
%       Ex: Numeric vectors
%
% PROPOSAL: Not require all "keys" to be unique for the same "key set"?
% PROPOSAL: Exclude empty strings (assertion).
%   PRO: Avoids empty string ambiguity.
%
% PROPOSAL: Return value as multiple return values (varying number).
%   NOTE: Requires values to be same sized 1D cell arrays.
%
% NOTE: Compare bicas.ga.normalize().

[keySetsTable, valuesTable] = convert_table(table);

% ASSERTION
for i = 1:numel(keySetsTable)
  assert(iscell(keySetsTable{i}), 'table(:, i) is not cell array.')
end

matchArray = zeros(size(keySetsTable));

nVarargin = numel(varargin);
if nVarargin == 1
  nonMatchError = false;
  nonMatchValue = varargin{1};
elseif nVarargin == 2
  nonMatchError = true;
  errorMsgId = varargin{1};
  errorMsg   = varargin{2};

  % NOTE: Empty errorMsgId ==> error() will not throw exception.
  assert(~isempty(errorMsgId), 'Empty errorMsgId')
  assert(ischar(errorMsgId))
  assert(ischar(errorMsg))
else
  error('Illegal number of arguments.')
end



% ASSERTIONS
% ASSERTION: Check for duplicate keys.
% NOTE: This condition requires the "keys" terms to be unique also within
% every set of keys.
combinedKeysList = [keySetsTable{:}];    % List of all keys in ALL RULES.
nCombinedKeys    = numel(combinedKeysList);
nUniqueKeys      = numel(unique(combinedKeysList));
assert(nCombinedKeys == nUniqueKeys, 'Illegal "table". Duplicated "keys" terms.')



%===========
% ALGORITHM
%===========
for i = 1:numel(keySetsTable)
  keySet = keySetsTable{i};

  % IMPLEMENTATION NOTE: ismember does not work as expected if keySet is
  % not a cell array, in particular a plain string.
  assert(iscell(keySet))

  matchArray(i) = ismember(key, keySet);
  if ismember(key, keySet)
    value = valuesTable{i};
    return
  end
end

% CASE: Did not find any match.
if nonMatchError
  error(errorMsgId, errorMsg)
else
  value = nonMatchValue;
end
end



function [keySetsTable, valuesTable] = convert_table(table)
irf.assert.sizes(table, [NaN, 2])
%assert(size(table, 2) == 2)

keySetsTable = table(:,1);
valuesTable  = table(:,2);

end
