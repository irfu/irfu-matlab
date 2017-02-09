% str = JSON_object_str(obj, settings)   Interprets a MATLAB variable as a JSON object and turns it into a string.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% NOTE: This is not a rigorous implementation based on a good understanding of JSON objects and therefore does not
% check for permitted characters.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% obj                : Recursively nested struct and cell arrays that can be interpreted as a JSON object.
%                      - Interprets MATLAB structure field names as "JSON parameter name strings".
%                      - Interprets MATLAB cell arrays as "JSON arrays/sets".
% settings           : Struct
%   .INDENT_SIZE     : Number of whitespace per indentation level.
%   .VALUE_POSITION  : The minimum number of characters between the beginning of a "name" and the beginning of the
%                      corresponding value. This setting can make the final string more readable.
% str                : Indented multi-line string with that is suitable for printing and human reading.
%                      NOTE: Uses line feed character for line breaks.
%
% NOTE: Since MATLAB structure field names are used for "JSON parameter name strings", the characters that can be used
% are likely more limited than what JSON permits.
%
function str = JSON_object_str(obj, settings)
%
% NOTE: Concerning the JSON syntax: "Whitespace is allowed and ignored around or between syntactic elements (values and
% punctuation, but not within a string value). Four specific characters are considered whitespace for this purpose:
% space, horizontal tab, line feed, and carriage return. JSON does not provide any syntax for comments."
% Source: https://en.wikipedia.org/wiki/JSON
%
% PROPOSAL: Permit arbitrary line break?

settings.LINE_BREAK = char(10);    % Should be same as sprintf('\n').

str = print_JSON_object_recursive(obj, 0, false, settings);
str = [str, settings.LINE_BREAK];

end



%###################################################################################################


% Recursive function that does the actual interpretation.
%
% indentFirstLine : True/false; Determines whether the first line should be indented.
%                   Can be useful for some layouts (placement of whitespace and line feeds).
function str = print_JSON_object_recursive(obj, indentationLevel, indentFirstLine, settings)

%LINE_BREAK = char(10);    % Should be same as sprintf('\n').
INDENT_0_STR = repmat(' ', 1, settings.INDENT_SIZE *  indentationLevel   );
INDENT_1_STR = repmat(' ', 1, settings.INDENT_SIZE * (indentationLevel+1));

str = '';


if iscell(obj)
    % CASE: Cell array - Interpret as a JSON array of (unnamed) JSON objects.
    
    if indentFirstLine
        str = [str, INDENT_0_STR];
    end
    str = [str, '[', settings.LINE_BREAK];
    
    for i = 1:length(obj)
        str = [str, print_JSON_object_recursive(obj{i}, indentationLevel+1, true, settings)];
        
        if i < length(obj)
            str = [str, ','];
        end
        str = [str, settings.LINE_BREAK];
    end
    str = [str, INDENT_0_STR, ']'];   % NOTE: No line break.
    
elseif isstruct(obj)
    % CASE: Struct - Interpret every field as value as (named) JSON object.
    
    if indentFirstLine
        str = [str, INDENT_0_STR];
    end
    str = [str, '{', settings.LINE_BREAK];
    
    names = fieldnames(obj);
    for i = 1:length(names)
        name = names{i};
        value = obj.(name);
        
        nameStr = sprintf('"%s":', name);
        str = [str, INDENT_1_STR, nameStr];
        
        if ischar(value)
            fillStr = repmat(' ', 1, settings.VALUE_POSITION-length(nameStr));
            str = [str, fillStr, sprintf('"%s"', value)];
        else
            % Alternative 1:
            % Left square brackets/braces begin on their own line.
            % Left & right brackets line up (same indentation). Is less (visually) compact.
            %str = [str, settings.LINE_BREAK, print_JSON_object_recursive(value, indentationLevel+1, true)];
            
            % Alternative 2:
            % Left square brackets/braces continue on the preceeding line.
            % Harder to visually match left & right brackets. More compact.
            str = [str, '  ', print_JSON_object_recursive(value, indentationLevel+1, false, settings)];
        end
        
        if i < length(names)
            str = [str, ','];
        end
        str = [str, settings.LINE_BREAK];
    end
    
    str = [str, INDENT_0_STR, '}'];   % NOTE: No line break.
    
else
    error('JSON_object_str:Assertion:IllegalArgument', 'Disallowed variable type. Neither structure nor cell array.')
end

end    % print_JSON_object_recursive
