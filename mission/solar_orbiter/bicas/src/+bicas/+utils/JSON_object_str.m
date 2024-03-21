%
% Interprets a MATLAB variable as a JSON object and turns it into a string for
% printing/writing to file.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% JsonObj
%       Recursively nested data structure that can be interpreted as a JSON
%       object.
%       - struct/containers.Map:
%           Field names/keys are interpreted as "JSON object names".
%           Values are interpreted recursively.
%       - cell array:
%           "JSON arrays/sets".
% indentSize
%       Number of whitespace per indentation level.
%
%
% RETURN VALUE
% ============
% str
%       Indented multi-line string that is suitable for printing and human
%       reading.
%       NOTE: Uses line feed character for line breaks.
%
%
% NOTES
% =====
% This is NOT a rigorous implementation based on a good understanding of JSON
% objects.
% It therefore
% ** does not check for permitted characters.
% ** likely does not support all the functionality of JSON:
%   ** JSON might permit more types of "leaves"/objects than supported by this
%      function.
% --
% Reasons to support:
%   struct         : Easier to hard-code. For backward compatibility.
%                    (Fieldnames are restricted. ==> object name strings are
%                    restricted.)
%   containers.Map : Can handle any object name string string.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-05-31
%
function str = JSON_object_str(JsonObj, indentSize)
%
% NOTE: Concerning the JSON syntax: "Whitespace is allowed and ignored
% around or between syntactic elements (values and punctuation, but not
% within a string value). Four specific characters are considered whitespace
% for this purpose: space, horizontal tab, line feed, and carriage return.
% JSON does not provide any syntax for comments."
% Source: https://en.wikipedia.org/wiki/JSON
%
% PROPOSAL: Automatic test.
%
% PROPOSAL: Permit arbitrary line break?

Settings.lineBreakStr = newline;    % Line feed
Settings.indentSize   = indentSize;

str = print_JSON_object_recursive(JsonObj, 0, false, Settings);
str = [str, Settings.lineBreakStr];
end



%###############################################################################



% Recursive function that does the actual interpretation.
%
% indentFirstLine
%       True/false; Determines whether the first line should be indented.
%       Can be useful for some layouts (placement of whitespace and line feeds).
function str = print_JSON_object_recursive(...
  JsonObj, indentationLevel, indentFirstLine, Settings)

%lineBreakStr = char(10);    % Should be same as sprintf('\n').
INDENT_0_STR = repmat(' ', 1, Settings.indentSize *  indentationLevel   );
INDENT_1_STR = repmat(' ', 1, Settings.indentSize * (indentationLevel+1));

str = '';



if iscell(JsonObj)
  % CASE: Cell array - Interpret as a JSON array of (unnamed) JSON objects.

  irf.assert.vector(JsonObj)

  if indentFirstLine
    str = [str, INDENT_0_STR];
  end
  str = [str, '[', Settings.lineBreakStr];

  for i = 1:length(JsonObj)
    % NOTE: RECURSIVE CALL
    str = [str, ...
      print_JSON_object_recursive(JsonObj{i}, ...
      indentationLevel+1, true, Settings)];

    if i < length(JsonObj)
      str = [str, ','];
    end
    str = [str, Settings.lineBreakStr];
  end
  str = [str, INDENT_0_STR, ']'];   % NOTE: No line break.

elseif isstruct(JsonObj) || isa(JsonObj, 'containers.Map')
  % CASE: Struct/containers.Map - Interpret every field/key as (named)
  %                               JSON object.

  [keysCa, valuesCa] = normalize_struct_Map_2_CA(JsonObj);

  if indentFirstLine
    str = [str, INDENT_0_STR];
  end
  str = [str, '{', Settings.lineBreakStr];

  maxLengthKey = max(cellfun(@length, keysCa));

  for i = 1:length(keysCa)
    key   = keysCa{i};
    value = valuesCa{i};

    keyStr = sprintf('"%s": ', key);   % NOTE: Four characters longer
    str = [str, INDENT_1_STR, keyStr];

    if ischar(value)
      % CASE: char value

      % NOTE: Hard-coded constant (aesthetics).
      fillStr = repmat(' ', 1, maxLengthKey+4 - length(keyStr));
      str     = [str, fillStr, sprintf('"%s"', value)];
    else
      % CASE: Non-char value ==> RECURSIVE CALL

      % Alternative 1:
      % Left square brackets/braces begin on their own line.
      % Left & right brackets line up (same indentation). Is less
      % (visually) compact.
      %                 str = [str, settings.lineBreakStr, ...
      %                     print_JSON_object_recursive(...
      %                         value, indentationLevel+1, true)];

      % Alternative 2:
      % Left square brackets/braces continue on the preceding line.
      % Harder to visually match left & right brackets. More compact.
      str = [str, '', print_JSON_object_recursive(...
        value, indentationLevel+1, false, Settings)];
    end

    if i < length(keysCa)
      % CASE: Not last key-value pair.
      str = [str, ','];
    end
    str = [str, Settings.lineBreakStr];
  end

  str = [str, INDENT_0_STR, '}'];   % NOTE: No line break.

else
  error('BICAS:Assertion:IllegalArgument', ...
    'Disallowed variable type. Neither structure nor cell array.')
end

end    % print_JSON_object_recursive



% Normalize struct and containers.Map to cell arrays (CA).
%
function [keysCa, valuesCa] = normalize_struct_Map_2_CA(JsonObj)
if isstruct(JsonObj)
  % CASE: struct

  keysCa   = fieldnames(JsonObj);
  valuesCa = cell(size(keysCa));
  for i = 1:numel(keysCa)
    valuesCa{i} = JsonObj.(keysCa{i});
  end

elseif isa(JsonObj, 'containers.Map')
  % CASE: containers.Map

  keysCa   = JsonObj.keys;
  valuesCa = cell(size(keysCa));
  for i = 1:numel(keysCa)
    assert(ischar(keysCa{i}))
    valuesCa{i} = JsonObj(keysCa{i});
  end

else
  error('Argument JsonObj is neither struct, nor containers.Map.')
end
end
