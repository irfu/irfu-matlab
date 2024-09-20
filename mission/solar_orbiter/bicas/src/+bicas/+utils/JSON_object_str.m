%
% Interprets a MATLAB variable as a JSON object and turns it into a string for
% printing/writing to file.
%
% NOTE: Can not set indentation length.
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
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-05-31
%
function str = JSON_object_str(JsonObj)
% PROPOSAL: Automatic test.

str = jsonencode(JsonObj, 'PrettyPrint', true);

% IMPLEMENTATION NOTE: String ends up in irf.str.add_prefix_on_every_row()
% which requires multi-row strings to always end with line feed.
str = [str, newline];    % Add line feed
end
