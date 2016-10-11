% select_structs   Given a cell array of structs, select a subset by matching (string) fields.
% 
% Generic function. Given a cell array of structs, select the subset of these where a specified field matches any of the
% specified string values.
%
%
% NOTE: Can be used multiple times successively to select for two matching struct fields, e.g. dataset_ID
% AND skeleton_version_str.
%
% ASSUMES: All structure field values (in "struct_list") for the chosen field ("field_name") are unique.
% ASSUMES: The structure field values are strings.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-19
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% struct_list        : Cell array of structs
% field_name         : Name of a struct field.
% selection_values   : Cell array of strings.
% subset_struct_list : Cell array of structs. That subset of structs in "struct_list" whose field name "field_name"
%                      equals any of the string values in "selection_values".
%
function subset_struct_list = select_structs(struct_list, field_name, selection_values)

% Convert to cell array of strings.
struct_list_values = {};
for i=1:length(struct_list)
    struct_list_values{i} = struct_list{i}.(field_name);
end

for i=1:length(selection_values)    
    j = find(strcmp(selection_values{i}, struct_list_values));
    if numel(j) ~= 1
        error('select_structs:Assertion:IllegalArgument', 'Can not find exactly one structure with field "%s" with value "%s".', field_name, selection_values{i})
    end
    subset_struct_list{i} = struct_list{j};
end

end
