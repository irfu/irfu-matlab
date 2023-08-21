%
% "Normalize" the fieldnames of a struct, i.e. assume that fieldnames can have
% multiple pre-defined variations (each actually used fieldname is one in a set)
% and convert these to the canonical fieldnames.
%
%
% NOTE: Could probably ~easily be generalized to struct arrays. The only(?) crux
% is that irf.ds.add_struct_to_struct() does not permit struct
% arrays.
% NOTE: Function can be used for always changing fieldname, since normList{i}{2}
% does not have to be a member of normList{i}{1}.
%
%
% ALGORITHM
% =========
% For every normalization rule:
% (1) Find those fieldnames which match the candidate list.
% (2) The first matching candidate field is renamed using the rule's designated
%     new fieldname.
% (3) Other matching candidates
%   (a) are ignored (fields are removed), or
%   (b) trigger error depending on policy.
%
%
% ARGUMENTS
% =========
% S1
%       Scalar struct.
% normList
%       List of "normalization rules".
%       {iRule}{1} = 1D cell array of strings. Candidate fieldnames which,
%                    when found in S, are (potentially) renamed.
%       {iRule}{2} = String. New fieldname. [] (numeric) == Remove field.
%       IMPLEMENTATION NOTE: Structure of normList is chosen to be suitable for hardcoding.
%       NOTE: The sets {iRule}{1} must not overlap.
% duplicatePolicy
%       String constant. Specificies assertion (what to permit).
%           '0-1 matches'
%           '0-inf matches'
%               Algorithm uses the first matching candidate and removes the rest.
%           '1 match,        'Assert one matching candidate'
%           '1-inf matches', 'Permit multiple matching candidates'
%               Algorithm uses the first matching candidate and removes the rest.
% --
%
%
% RETURN VALUES
% =============
% S2
%       Modified argument struct.
% changeList
%       Struct array of fieldname changes. Summary of fieldname changes that
%       actually took place.
%           .oldFieldname
%               Old fieldname.
%           .newFieldname
%               New fieldname. Is never identical to the corresponding
%               .oldFieldname.
%               [] (numeric) == Removed field.
%           .ignoredCandidateFieldnames
%               Cell array of strings.
%       NOTE: This is useful for the caller to vary its behaviour depending on
%             changes.
%             Ex: Choose between doing (a) nothing, (b) warning, (c) error, and
%                 behaving differently depending on fieldname (e.g. error for
%                 one change, but warning for another).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-12-02.
%
function [S2, fnChangeList] = normalize_struct_fieldnames(...
  S1, normList, duplicatePolicy)
%
% NOTE: May eventually be implemented using irf.utils.translate_strings.
% NOTE: See BOGIQ for                       irf.utils.translate_strings.
%
% NOTE: Can not be used to directly handle zVariable names in dataobj, since zVariable names are stored in multiple
% locations.
%   Ex: irfu-matlab code uses zVarnames in other places than DataObj.data fieldnames.
%       Ex: getfillval, findva (used by getfillval).
%
% TODO-DEC: How handle fieldnames which are not found?
%   PROPOSAL: Better to "translate" (change concept) when fieldnames are found, and ignore otherwise.
%       PRO/NOTE: Can already "translate" fieldnames (change names arbitrarily) since normList{i}{2} does not have
%                 to be a member of normList{i}{1}.
%   PROPOSAL: Policy arguments
%       Alt 1: 0-1   matches per rule
%       Alt 2: 0-inf matches per rule
%       Alt 3: 1     match per rule
%       Alt 4: 1-inf matches per rule
%           NOTE: Permitting zero matches is NOT the same as including the new fieldname among the candidates 8and
%           requiring <=1 matches), since this requires the field to exist.
%       PROPOSAL: Use irf.utils.interpret_settings_args.
%
% PROBLEM: Need to rethink how changes are represented. ignoredCandidateFieldnames now contains fields which may
% have been removed.
%   NOTE: This would also be useful if re-writing to normalize set of (unique) strings.
%   NEED: Enable customized handling depending on execution.
%       Ex: Ignore, warning, error depending on whether certain rule changes/removes a fieldname.
%       Ex: Customized log messages depending on what happened.
%       Ex: bicas.handle_struct_name_change
%   NEED: Way of reversing changes(?)
%       Ex: When modifying dataset: Normalize struct fieldnames, modify struct, de-normalize struct fieldnames.
%   PROPOSAL: Represent changes as what happens to each original field.
%       (1) Kept
%       (2) Renamed: new fieldname
%       (3) Removed
%       PROPOSAL: Table: Changes(iField)
%           .oldFieldname : Note: Special value could represent adding/creating (did not pre-exist).
%           .newFieldname : Special value represents removal
%       PROPOSAL: Document all correlations between before-after. Similar to permutations, sort().
%           i2 = i12(i1), i1 = i21(i2)
%           NOTE: Still needs to represent deletions.
%           PRO: Good for denormalization.
%   PROPOSAL: Represent changes as how every rule was acted upon.
%       (1) Which candidate was used.
%       (2) Which candidates were ignored.
%
% PROBLEM: Concept of changing fieldnames of actual struct is bad when struct is so large that it may cause memory
%          problems. May want to avoid modifying struct in order to help MATLAB's code optimization (prevent
%          temporary copies).

fnList       = fieldnames(S1);
iChange      = 0;
fnChangeList = irf.ds.empty_struct(...
  [0,1], 'oldFieldname', 'newFieldname', 'ignoredCandidateFieldnames');
% IMPLEMENTATION NOTE: Initializing with empty struct mean that calling code
% does not need a special case for e.g. {fnChangeList.oldFieldname}.

% NR = Normalization Rule
S2 = irf.ds.empty_struct(size(S1));
% Iterate over items in NORMALIZATION LIST (not over fieldnames in S).
for iNr = 1:numel(normList)
  normRule = normList{iNr};

  % ASSERTION: Check argument format.
  assert(numel(normRule) == 2, ...
    'Not exactly two elements in argument normList{%i}.', iNr)

  fnCandidates = normRule{1};
  fnNew        = normRule{2};

  % ASSERTIONS
  assert(~isempty(fnCandidates), ...
    'Argument normList{%i}{1} is empty.', iNr)



  %=========================================================
  % Find (1) candidate to use, and (2) candidates to ignore
  %=========================================================
  % Find indices into fnCandidates (not fnList): All candidates which
  % match a field name (one, several).
  iFn = find(ismember(fnCandidates, fnList));
  ignoredCandidateFieldnames = cell(1,0);
  switch(duplicatePolicy)

    case '0-1 matches'
      if isempty(iFn)
        continue
      else
        % ASSERTION
        if ~isscalar(iFn)
          policy_error(...
            'Did not find zero or one field name that matches any of %s.', ...
            iNr, fnCandidates)
        end
      end

    case '0-inf matches'
      if isempty(iFn)
        continue
      end

      % Keep the first matching candidate, list the ignored ones.
      ignoredCandidateFieldnames = fnCandidates(iFn(2:end));
      iFn = iFn(1);

    case {'1 match', 'Assert one matching candidate'}
      % ASSERTION
      if ~isscalar(iFn)
        policy_error(...
          'Did not find exactly one field name that matches any of %s', ...
          iNr, fnCandidates)
      end

    case {'1-inf matches', 'Permit multiple matching candidates'}
      % ASSERTION
      if ~(numel(iFn) >= 1)
        policy_warning(...
          'Did not find at least on field name that matches any of %s.', ...
          iNr, fnCandidates)
      end

      % Keep the first matching candidate, list the ignored ones.
      ignoredCandidateFieldnames = fnCandidates(iFn(2:end));
      iFn = iFn(1);

    otherwise
      error('Ilegal argument duplicatePolicy="%s".', duplicatePolicy)
  end
  fnOld = fnCandidates{iFn};

  %==================================================
  % "Execute" and "log" the resulting change, if any
  %==================================================
  if ischar(fnNew) && ~isempty(fnNew)
    S2.(fnNew) = S1.(fnOld);     % COPY field
  elseif isnumeric(fnNew) && isempty(fnNew)
    % Do nothing. Ignore/remove old field.
  else
    error('Illegal argument normList{%i}{2}', iNr)
  end

  % REMOVE fields
  S1 = rmfield(S1, [{fnOld}, ignoredCandidateFieldnames(:)']);

  if ~strcmp(fnOld, fnNew) || ~isempty(ignoredCandidateFieldnames)
    iChange = iChange + 1;
    fnChangeList(iChange, 1).ignoredCandidateFieldnames = ignoredCandidateFieldnames;
    fnChangeList(iChange, 1).oldFieldname               = fnOld;
    fnChangeList(iChange, 1).newFieldname               = fnNew;
  end
end

% Move fields which were not affected by normalization rules.
% NOTE: add_struct_to_struct() does work with arrays, but only if it does
% not need to inspect the field values.
S2 = irf.ds.add_struct_to_struct(S2, S1);
end



% msg : String containing one occurrence of "%s" for the candidate fieldnames.
function policy_error(msg, iNr, fnCandidates)
normListStr = sprintf('normList{%i}{1}={"%s"}.', ...
  iNr, strjoin(fnCandidates, '", "'));
msg = sprintf(msg, normListStr);
error(msg)
end

