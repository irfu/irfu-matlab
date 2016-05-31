%
% Interprets structure as a JSON object and returns it as string that could be printed.
%
% - Interprets MATLAB stucture field names as JSON parameter name strings.
%   NOTE: Does in principle limit the characters that can be used for such.
% - Interprets MATLAB cell arrays as JSON arrays.
%
function str = JSON_object_str(obj)

INDENT_SIZE = 4;
FILL_STR_MAX_LENGTH = 13;
LF = sprintf('\n');


str = print_JSON_object_recursive(obj, 0, true);
str = [str, LF];

%==========================================================================
    % indent_first_line = True/false; Determines whether the first line should be indented.
    % Can be useful for some layouts (placement of whitespace and line feeds).
    function str = print_JSON_object_recursive(obj, indent_level, indent_first_line)
        
        str = '';
        
        indent0_str = repmat(' ', 1, INDENT_SIZE* indent_level);
        indent1_str = repmat(' ', 1, INDENT_SIZE*(indent_level+1));
        
        if iscell(obj)
            % CASE: Interpret cell value as a JSON array of nested JSON objects.
            
            if indent_first_line
                str = [str, indent0_str];
            end
            str = [str, '[', LF];
            
            for i = 1:length(obj)
                str = [str, print_JSON_object_recursive(obj{i}, indent_level+1, true)];
                
                if i < length(obj)
                    str = [str, ','];
                end
                str = [str, LF];
            end
            str = [str, sprintf('%s]', indent0_str)];

        else
            % CASE: Interpret the value as a non-array JSON object.
            
            if indent_first_line
                str = [str, indent0_str];
            end
            str = [str, '{', LF];
            
            names = fieldnames(obj);
            for i = 1:length(names)
                name = names{i};
                value = obj.(name);
                
                str = [str, indent1_str, sprintf('"%s":', name)];
                
                if ischar(value)
                    fill_str = repmat(' ', 1, FILL_STR_MAX_LENGTH-length(name));
                    str = [str, fill_str, sprintf('"%s"', value)];
                else
                    % Alternative 1: Left square brackets/braces begin on their own line.
                    % Opening & ending brackets line up (same indentation). Is less compact.
                    %str = [str, LF, print_JSON_object_recursive(value, indent_level+1, true)];
                    
                    % Alternative 2: Left square brackets/braces continue on the preceeding line.
                    % Harder to match opening/closing brackets. More compact.
                    str = [str, '  ', print_JSON_object_recursive(value, indent_level+1, false)];
                end
                
                if i < length(names)
                    str = [str, ','];                    
                end
                str = [str, LF];
            end
            
            str = [str, indent0_str, '}'];
            
        end
    end    % print_JSON_object_recursive
end
