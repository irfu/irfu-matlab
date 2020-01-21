%
% "Normalize" the fieldnames of a struct, i.e. assume that fieldnames can have multiple pre-defined variations (each
% actually used fieldname is one in a set) and convert these to the canonical fieldnames.
%
%
% ARGUMENTS
% =========
% S        : Arbitrary struct
% normList : "Normalization list"
%            {i}{1} = Cell array of fieldnames which, when found in S, are renamed.
%            {i}{2} = (a string) New fieldname.
%            IMPLEMENTATION NOTE: Structure of normList is chosen to be suitable for hardcoding.
%            NOTE: The sets {i}{1} must not overlap.
% --
% Can also be used for always changing fieldname, since normList{i}{2} does not have to be a member of normList{i}{1}.
%
%
% RETURN VALUES
% =============
% S          : Modified argument struct.
% changeList : Struct array of fieldname changes. Summary of fieldname changes that actually took place.
%                 .oldFieldname = Old fieldname
%                 .newFieldname = New fieldname. Is never identical to the corresponding .oldFieldname.
%                 NOTE: This is useful for the caller vary its behaviour depending on changes.
%                   Ex: Choose between doing (a) nothing, (b) warning, (c) error, and beahving differently depending on
%                   fieldname (e.g. error for one change, but warning for another).
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2019-12-02.
%
function [S, fnChangeList] = normalize_struct_fieldnames(S, normList)
    % TODO-DECISION: How handle fieldnames which are not found?
    %   PROPOSAL: Better to "translate" (change concept) when fieldnames are found, and ignore otherwise.
    %       PRO/NOTE: Can already "translate" fieldnames (change names arbitrarily) since normList{i}{2} does not have to be a member of normList{i}{1}.
    %   PROPOSAL: Policy argument.
    % PROPOSAL: Test code.
    %
    % PROPOSAL: Return ~flags for what has been changed.
    %   PRO: Caller can print warnings based on it.

    fnList       = fieldnames(S);
    iChange      = 0;
    fnChangeList = struct('oldFieldname', {}, 'newFieldname', {});
    % IMPLEMENTATION NOTE: Initializing with empty struct meand that calling code does not need a special case for e.g.
    % {fnChangeList.oldFieldname}.

    % NL = Normalization List
    for iNl = 1:numel(normList)    % Iterate over items in NORMALIZATION LIST (not over fieldnames in S).
        % ASSERTION: Check argument format.
        assert(numel(normList{iNl}) == 2, 'Not exactly two elements in normList{%i}.', iNl)
                
        fnOldAlternatives = normList{iNl}{1};
        fnNew             = normList{iNl}{2};
        
        % ASSERTIONS
        assert(~isempty(fnOldAlternatives), 'normList{%i}{1} is empty.', iNl)
        EJ_library.assert.castring(fnNew)
        
        % Find indices into fnList: Fieldnames which match normalized field name iNl (should be exactly one).
        iFn = find(ismember(fnList, fnOldAlternatives));
        % ASSERTION: There is exactly one fieldname that matches the current canonical fieldname (item in normList).
        assert(isscalar(iFn), 'Did not find exactly one field name that matches any of normList{%i}{1}={"%s"}.', iNl, strjoin(fnOldAlternatives, '", "'))

        fnOld = fnList{iFn};

        if ~strcmp(fnOld, fnNew)
            iChange = iChange + 1;
            fnChangeList(iChange).oldFieldname = fnOld;
            fnChangeList(iChange).newFieldname = fnNew;
            S.(fnNew)                          = S.(fnOld);
            S                                  = rmfield(S, fnOld);
        end
    end
end
