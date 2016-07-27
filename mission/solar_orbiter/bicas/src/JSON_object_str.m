% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Interprets a MATLAB variable consisting of recursively nested structures and cell arrays as a JSON
% object and returns it as an indented multi-line string that is suitable for printing and human
% reading.
%
% - Interprets MATLAB structure field names as JSON parameter name strings.
%   NOTE: Does in principle limit the characters that can be used for such.
% - Interprets MATLAB cell arrays as JSON arrays.
%
% NOTE: This is not a rigorous implementation and therefore does not check for permitted characters.
%
function str = JSON_object_str(obj, settings)
%
% NOTE: "Whitespace is allowed and ignored around or between syntactic elements (values and
% punctuation, but not within a string value). Four specific characters are considered whitespace
% for this purpose: space, horizontal tab, line feed, and carriage return. JSON does not provide
% any syntax for comments."
% https://en.wikipedia.org/wiki/JSON



LINE_BREAK = sprintf('\n');    % Put in constants structure?!!

str = print_JSON_object_recursive(obj, 0, false, settings);
str = [str, LINE_BREAK];

end



%###################################################################################################



% indent_first_line = True/false; Determines whether the first line should be indented.
% Can be useful for some layouts (placement of whitespace and line feeds).
function str = print_JSON_object_recursive(obj, indent_level, indent_first_line, settings)

LINE_BREAK = sprintf('\n');    % Put in constants structure?!!
INDENT_0_STR = repmat(' ', 1, settings.indent_size *  indent_level   );
INDENT_1_STR = repmat(' ', 1, settings.indent_size * (indent_level+1));

str = '';


if iscell(obj)
    % CASE: Cell array - Interpret as a JSON array of (unnamed) JSON objects.
    
    if indent_first_line
        str = [str, INDENT_0_STR];
    end
    str = [str, '[', LINE_BREAK];
    
    for i = 1:length(obj)
        str = [str, print_JSON_object_recursive(obj{i}, indent_level+1, true, settings)];
        
        if i < length(obj)
            str = [str, ','];
        end
        str = [str, LINE_BREAK];
    end
    str = [str, INDENT_0_STR, ']'];   % NOTE: No line break.
    
elseif isstruct(obj)
    % CASE: Struct - Interpret every field as value as (named) JSON object.
    
    if indent_first_line
        str = [str, INDENT_0_STR];
    end
    str = [str, '{', LINE_BREAK];
    
    names = fieldnames(obj);
    for i = 1:length(names)
        name = names{i};
        value = obj.(name);
        
        str = [str, INDENT_1_STR, sprintf('"%s":', name)];
        
        if ischar(value)
            fill_str = repmat(' ', 1, settings.fill_str_max_length-length(name));
            str = [str, fill_str, sprintf('"%s"', value)];
        else
            % Alternative 1:
            % Left square brackets/braces begin on their own line.
            % Left & right brackets line up (same indentation). Is less compact.
            %str = [str, LINE_BREAK, print_JSON_object_recursive(value, indent_level+1, true)];
            
            % Alternative 2:
            % Left square brackets/braces continue on the preceeding line.
            % Harder to visually match left & right brackets. More compact.
            str = [str, '  ', print_JSON_object_recursive(value, indent_level+1, false, settings)];
        end
        
        if i < length(names)
            str = [str, ','];
        end
        str = [str, LINE_BREAK];
    end
    
    str = [str, INDENT_0_STR, '}'];   % NOTE: No line break.
    
else
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Disallowed variable type. Neither structure nor cell array.')
end

end    % print_JSON_object_recursive
