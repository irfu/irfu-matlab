% Interprets structure as a JSON object and prints it to stdout as such.
%
% - Interprets MATLAB field names as JSON parameter name strings.
%   NOTE: Does in principle limit the characters that can be used for such.
% - Interprets MATLAB cell arrays as JSON arrays.
%
function print_JSON_object(obj)
print_JSON_object_recursive(obj, 0);
fprintf('\n');

%==========================================================================
    function print_JSON_object_recursive(obj, indent_level)
        INDENT_SIZE = 4;
        FILL_STR_MAX_LENGTH = 12;
        
        indent0_str = repmat(' ', 1, INDENT_SIZE* indent_level);
        indent1_str = repmat(' ', 1, INDENT_SIZE*(indent_level+1));
        
        if iscell(obj)
            % CASE: Interpret cell value as JSON array of nested JSON objects.
            
            fprintf('[\n%s', indent1_str)
            for i = 1:length(obj)
                print_JSON_object_recursive(obj{i}, indent_level+1)
                
                if i < length(obj)
                    fprintf(',');
                end
                fprintf('\n');
            end
            fprintf('%s]', indent0_str)
            
        else
            % CASE: Interpret the value as a non-array JSON object.
            
            % No indentation since it is assumed that there is a preceeding name string
            % except for the outermost object which should not be indented anyway.
            fprintf('{\n');
            
            names = fieldnames(obj);
            for i = 1:length(names)
                name = names{i};
                value = obj.(name);
                
                fill_str = repmat(' ', 1, FILL_STR_MAX_LENGTH-length(name));
                fprintf('%s"%s": %s', indent1_str, name, fill_str)
                
                if ischar(value)
                    fprintf('"%s"', value)
                else
                    print_JSON_object_recursive(value, indent_level+1);
                end
                
                if i < length(names)
                    fprintf(',');
                end
                fprintf('\n');
            end
            
            fprintf('%s}', indent0_str);
            
        end
    end    % print_JSON_object_recursive
end
