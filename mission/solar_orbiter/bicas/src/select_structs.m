% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-19
%
% Generic function.
% Given a cell array of structs, select a subset of these structs based upon the (string) values of
% one of their fields.
%
% NOTE: Can be used multiple times successively to select using two struct fields, e.g. dataset_ID
% AND dataset_version_str.
%
% ASSUMES: All structure field values (in struct_list) for the chosen field are unique.
% ASSUMES: The structure field values are strings.
%
% struct_list = cell array of structs
% field_name = string
% field_valies = cell array of strings
% return_stuct_list = cell array of structs
%
function return_struct_list = select_structs(struct_list, field_name, selection_values)
global ERROR_CODES

% Convert to cell array of strings.
struct_list_values = {};
for i=1:length(struct_list)
    struct_list_values{i} = struct_list{i}.(field_name);
end

for i=1:length(selection_values)    
    j = find(strcmp(selection_values{i}, struct_list_values));
    if numel(j) ~= 1
        errorp(ERROR_CODES.ASSERTION_ERROR, 'Pure code bug. Can not find exactly one structure for field value "%s".', selection_values{i})
    end
    return_struct_list{i} = struct_list{j};
end

end

