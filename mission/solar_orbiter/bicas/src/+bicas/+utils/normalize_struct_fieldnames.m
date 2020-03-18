%
% "Normalize" the fieldnames of a struct, i.e. assume that fieldnames can have multiple pre-defined variations (each
% actually used fieldname is one in a set) and convert these to the canonical fieldnames.
%
%
% NOTE: Could probably ~easily be generalized to struct arrays. The only(?) crux is that
% EJ_library.utils.add_struct_to_struct does not permit struct arrays.
% NOTE: Function can be used for always changing fieldname, since normList{i}{2} does not have to be a member of
% normList{i}{1}.
%
%
% ALGORITHM
% =========
% For every normalization rule:
% (1) Find those fieldnames which match the candidate list.
% (2) The first matching candidate field is renamed using the rule's designated new fieldname.
% (3) Other matching candidates are ignored (removed), or trigger error depending on policy.
%
%
% ARGUMENTS
% =========
% S1       : Scalar struct.
% normList : List of "normalization rules".
%            {iRule}{1} = 1D cell array of strings. Candidate fieldnames which, when found in S, are (potentially) renamed.
%            {iRule}{2} = String. New fieldname.
%            IMPLEMENTATION NOTE: Structure of normList is chosen to be suitable for hardcoding.
%            NOTE: The sets {iRule}{1} must not overlap.
% duplicatePolicy : String constant. 
%           'Assert one matching candidate'
%           'Permit multiple matching candidates'. Algorithm uses the first matching candidate and removes the rest.
% --
%
%
% RETURN VALUES
% =============
% S2         : Modified argument struct.
% changeList : Struct array of fieldname changes. Summary of fieldname changes that actually took place.
%   .oldFieldname               : Old fieldname.
%   .newFieldname               : New fieldname. Is never identical to the corresponding .oldFieldname.
%   .ignoredCandidateFieldnames : Cell array of strings.
%                 NOTE: This is useful for the caller vary its behaviour depending on changes.
%                 Ex: Choose between doing (a) nothing, (b) warning, (c) error, and beahving differently depending on
%                     fieldname (e.g. error for one change, but warning for another).
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2019-12-02.
%
function [S2, fnChangeList] = normalize_struct_fieldnames(S1, normList, duplicatePolicy)
    %
    % TODO-DECISION: How handle fieldnames which are not found?
    %   PROPOSAL: Better to "translate" (change concept) when fieldnames are found, and ignore otherwise.
    %       PRO/NOTE: Can already "translate" fieldnames (change names arbitrarily) since normList{i}{2} does not have
    %                 to be a member of normList{i}{1}.
    %   PROPOSAL: Policy argument ~unusedRulePolicy, ~unusedNormalizationRulePolicy.
    %       Alt 1: Ignore.
    %       Alt 2: Assert every item in normalization list is used.
    %       PROPOSAL: Use EJ_library.utils.interpret_settings_args.
    %
    % PROPOSAL: Rewrite into function for normalizing set of strings (corresponding to set of fieldnames in struct now).
    %   This function can then be implemented using the new function, assuming it returns a good enough change list.
    %   PRO: Such a function would be somewhat easier to automatically test.
    %   PRO: Such a function would potentially have a wider range of uses.
    %       CON: Can't think of any examples?
    %   NOTE: Cf EJ_library.utils.translate.

    fnList       = fieldnames(S1);
    iChange      = 0;
    fnChangeList = EJ_library.utils.empty_struct([0,1], 'oldFieldname', 'newFieldname', 'ignoredCandidateFieldnames');
    % IMPLEMENTATION NOTE: Initializing with empty struct mean that calling code does not need a special case for e.g.
    % {fnChangeList.oldFieldname}.

    % NL = Normalization List
    S2 = EJ_library.utils.empty_struct(size(S1));
    for iNl = 1:numel(normList)    % Iterate over items in NORMALIZATION LIST (not over fieldnames in S).
        % ASSERTION: Check argument format.
        assert(numel(normList{iNl}) == 2, 'Not exactly two elements in argument normList{%i}.', iNl)
                
        fnCandidates = normList{iNl}{1};
        fnNew        = normList{iNl}{2};
        
        % ASSERTIONS
        assert(~isempty(fnCandidates), 'Argument normList{%i}{1} is empty.', iNl)
        EJ_library.assert.castring(fnNew)
        
        
        
        %=========================================================
        % Find (1) candidate to use, and (2) candidates to ignore
        %=========================================================
        % Find indices into fnCandidates (not fnList): All candidates which match a field name (one, several).
        iFn = find(ismember(fnCandidates, fnList));
        ignoredCandidateFieldnames = cell(1,0);
        switch(duplicatePolicy)

            case 'Assert one matching candidate'
                % ASSERTION: There is exactly one fieldname that matches the current canonical fieldname (item in normList).
                assert(isscalar(iFn), ...
                    'Did not find exactly one field name that matches any of normList{%i}{1}={"%s"}.', ...
                    iNl, strjoin(fnCandidates, '", "'))

            case 'Permit multiple matching candidates'
                % Keep the first matching candidate, list the ignored ones.
                ignoredCandidateFieldnames = fnCandidates(iFn(2:end));
                iFn = iFn(1);

            otherwise
                error('Ilegal argument duplicatePolicy="%s".', duplicatePolicy)
        end
        fnOld = fnCandidates{iFn};

        %====================================================
        % "Implement" and "log" the resulting change, if any
        %====================================================
        S2.(fnNew) = S1.(fnOld);
        S1 = rmfield(S1, [{fnOld}, ignoredCandidateFieldnames(:)']);
        if ~strcmp(fnOld, fnNew) || ~isempty(ignoredCandidateFieldnames)
            iChange = iChange + 1;
            fnChangeList(iChange, 1).ignoredCandidateFieldnames = ignoredCandidateFieldnames;
            fnChangeList(iChange, 1).oldFieldname               = fnOld;
            fnChangeList(iChange, 1).newFieldname               = fnNew;
        end
    end
    
    % Move fields which were not affected by normalization rules.
    S2 = EJ_library.utils.add_struct_to_struct(S2, S1);
end
