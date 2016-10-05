% Assert that a cell 1D array of strings only contains unique strings.
%
% This is useful for strings which serve as "identifiers", that can be "pointed to" and have to be unique.
%
function assert_strings_unique(string_list)
%
% PROPOSAL: Have caller throw error instead. Return list of unique strings.
% PROPOSAL: Generalize to arbitrary data using ~isequal, isequaln.
%
[unique_strings, i_l, i_u] = unique(string_list);    % iu=indices to unique_strings; il=indices to string_list

if length(unique_strings) == length(string_list)
    %return
    ;
else
    % Figure out which strings occur more than one time. Surprisingly difficult for such a simple
    % task.
    % "TEST CODE" for playing around on the command line.
    % l={'a', 'qwe', 'asd', 'qwe', 'qwe', 'a'} ; [u, i_l, i_u]=unique(l); i_us = sort(i_u); i_i_us = find(~diff(i_us)); unique(u(i_us(i_i_us)))
    % ----------------------------------------------------------------------
    % i_u : Each unique component value represents a unique string.
    i_us = sort(i_u);                % s=Sorted "strings"
    i_i_us = find(~diff(i_us));      % Indices into i_us for component values ("strings") in i_us occurring twice in a row.
    nonuniques_list = unique(unique_strings(i_us(i_i_us)));

    
    disp_list = '';
    for i = 1:length(nonuniques_list)
        disp_list = [disp_list, sprintf('   "%s"\n', nonuniques_list{i})];
    end    
    
    error('assert_strings_unique:Assertion', ...
        ['Strings in list of strings are not unique as expected. The following quoted strings are not unique:\n', disp_list])
end

end
