% Interprets a MATLAB variable as a JSON object and turns it into a string for printing/writing to file.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% jsonObj            : Recursively nested struct and cell arrays that can be interpreted as a JSON object.
%                      - Interprets MATLAB structure field names as "JSON parameter name strings".
%                      - Interprets MATLAB cell arrays as "JSON arrays/sets".
% indentSize         : Number of whitespace per indentation level.
% valuePosition      : The minimum number of characters between the beginning of a "name" and the beginning of the
%                      corresponding value. This setting can make the final string more readable.
% str                : Indented multi-line string that is suitable for printing and human reading.
%                      NOTE: Uses line feed character for line breaks.
%
%
% NOTE: This is NOT a rigorous implementation based on a good understanding of JSON objects and therefore does not
% check for permitted characters.
%
% NOTE: Since MATLAB structure field names are used for "JSON parameter name strings", the characters that can be used
% are likely more limited than what JSON permits.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
function str = JSON_object_str(JsonObj, indentSize)
%
% NOTE: Concerning the JSON syntax: "Whitespace is allowed and ignored around or between syntactic elements (values and
% punctuation, but not within a string value). Four specific characters are considered whitespace for this purpose:
% space, horizontal tab, line feed, and carriage return. JSON does not provide any syntax for comments."
% Source: https://en.wikipedia.org/wiki/JSON
%
% PROPOSAL: Permit arbitrary line break?

Settings.LINE_BREAK     = char(10);    % Line feed
Settings.INDENT_SIZE    = indentSize;

str = print_JSON_object_recursive(JsonObj, 0, false, Settings);
str = [str, Settings.LINE_BREAK];

end



%###################################################################################################


% Recursive function that does the actual interpretation.
%
% indentFirstLine : True/false; Determines whether the first line should be indented.
%                   Can be useful for some layouts (placement of whitespace and line feeds).
function str = print_JSON_object_recursive(JsonObj, indentationLevel, indentFirstLine, Settings)

%LINE_BREAK = char(10);    % Should be same as sprintf('\n').
INDENT_0_STR = repmat(' ', 1, Settings.INDENT_SIZE *  indentationLevel   );
INDENT_1_STR = repmat(' ', 1, Settings.INDENT_SIZE * (indentationLevel+1));

str = '';



if iscell(JsonObj)
    % CASE: Cell array - Interpret as a JSON array of (unnamed) JSON objects.
    
    if indentFirstLine
        str = [str, INDENT_0_STR];
    end
    str = [str, '[', Settings.LINE_BREAK];
    
    for i = 1:length(JsonObj)
        str = [str, print_JSON_object_recursive(JsonObj{i}, indentationLevel+1, true, Settings)];    % NOTE: RECURSIVE CALL
        
        if i < length(JsonObj)
            str = [str, ','];
        end
        str = [str, Settings.LINE_BREAK];
    end
    str = [str, INDENT_0_STR, ']'];   % NOTE: No line break.
    
elseif isstruct(JsonObj)
    % CASE: Struct - Interpret every field as value as (named) JSON object.
    
    if indentFirstLine
        str = [str, INDENT_0_STR];
    end
    str = [str, '{', Settings.LINE_BREAK];
    
    fnList = fieldnames(JsonObj);
    maxLengthFn = max(cellfun(@length, fnList));
    
    for i = 1:length(fnList)
        name = fnList{i};
        value = JsonObj.(name);
        
        nameStr = sprintf('"%s": ', name);   % NOTE: Four characters longer
        str = [str, INDENT_1_STR, nameStr];
        
        if ischar(value)
            % CASE: char value
            fillStr = repmat(' ', 1, maxLengthFn+4 - length(nameStr));    % NOTE: Hard-coded constant (estethics).
            str = [str, fillStr, sprintf('"%s"', value)];
        else
            % CASE: Non-char value ==> recurse
            
            % Alternative 1:
            % Left square brackets/braces begin on their own line.
            % Left & right brackets line up (same indentation). Is less (visually) compact.
            %str = [str, settings.LINE_BREAK, print_JSON_object_recursive(value, indentationLevel+1, true)];
            
            % Alternative 2:
            % Left square brackets/braces continue on the preceeding line.
            % Harder to visually match left & right brackets. More compact.
            str = [str, '', print_JSON_object_recursive(value, indentationLevel+1, false, Settings)];
        end
        
        if i < length(fnList)
            str = [str, ','];
        end
        str = [str, Settings.LINE_BREAK];
    end
    
    str = [str, INDENT_0_STR, '}'];   % NOTE: No line break.
    
else
    error('BICAS:JSON_object_str:Assertion:IllegalArgument', 'Disallowed variable type. Neither structure nor cell array.')
end

end    % print_JSON_object_recursive
