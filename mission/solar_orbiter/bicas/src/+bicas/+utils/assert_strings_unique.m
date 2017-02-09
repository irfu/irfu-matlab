% assert_strings_unique(stringList)   Assert that a cell array of strings only contains unique strings.
%
% This is useful for checking lists of strings which serve as "identifiers", that can be "pointed to" and that have to
% be unique.
%
function assert_strings_unique(stringList)
%
% PROPOSAL: Have caller throw error instead. Return list of unique strings.
% PROPOSAL: Generalize to arbitrary data using ~isequal, isequaln.

[uniqueStringsList, i_l, i_u] = unique(stringList);    % i_u=indices to unique_strings; i_l=indices to string_list

if numel(uniqueStringsList) == numel(stringList)
    ;   % Do nothing.
else
    % Figure out which strings that occur more than one time. This is surprisingly difficult for such a simple task.
    % "TEST CODE" for playing around on the command line.
    % l={'a', 'qwe', 'asd', 'qwe', 'qwe', 'a'} ; [u, i_l, i_u]=unique(l); i_us = sort(i_u); i_i_us = find(~diff(i_us)); unique(u(i_us(i_i_us)))
    % ----------------------------------------------------------------------
    % i_u : Each unique component value represents a unique string.
    i_us = sort(i_u);                % s=Sorted "strings"
    i_i_us = find(~diff(i_us));      % Indices into i_us for component values ("strings") in i_us occurring twice in a row.
    nonuniquesList = unique(uniqueStringsList(i_us(i_i_us)));

    
    dispList = '';
    for i = 1:numel(nonuniquesList)
        dispList = [dispList, sprintf('   "%s"\n', nonuniquesList{i})];
    end    
    
    error('assert_strings_unique:Assertion', ...
        ['Strings in list of strings are not unique as expected. The following quoted strings are not unique:\n', dispList])
end

end
