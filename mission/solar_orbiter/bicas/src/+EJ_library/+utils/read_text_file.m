%
% Reads file as text file and stores lines in cell array.
%
%
% RETURN VALUE
% ============
% rowsList : Column cell array of strings (one per row). All characters in the line breaks have been removed so
%            that it works for both Windows-style and Unix-style text files (assuming linebreakRegexp has been set
%            correctly).
%            NOTE: There will be a row after the last line break. The caller has to decide whether this is a legitimate
%            row.
% linebreakRegexp : Regular expression used as line break.
%
%
% IMPLEMENTATION NOTE
% ===================
% Old implementation using textscan merged consequtive linebreaks (CR+LF) into one (i.e. bug). Therefore not using.
%
function [rowsList] = read_text_file(filePath, linebreakRegexp)
% PROPOSAL: Use s=textscan(fileId, '%s', 'delimiter', '\n').
% PROPOSAL: Read entire file as a string, then split it using a specified line-break string.
%   PROPOSAL: Return both string and row string list.
%       CON: Unnecessary to return row string list. Trivial to produce with str_split.
%           CON: Only if what is line-break is unambiguous.
%       CON: Row string list is line break-dependent.
%           CON: If want automatic detection of type of line break, then rowsList is the proper format, not the
%                multi-row string.
%       CON: Row string list is dependent on how to interpret chars after the last line break.
%
% PROPOSAL: Split into functions
%   (1) read entire file as sequence of bytes.
%   (2) convert sequence of bytes into rows
%   PRO: Can have test code for converting sequence into rows.
%   PRO: Can have different algorithms (functions) for different conversion.

    fileId = fopen(filePath);
    
    % ASSERTION
    if fileId == -1
        error('read_text_file:CanNotOpenFile', 'Can not open file: "%s"', filePath)
    end
    
% 
%     % IMPLEMENTATION NOTE: "delimiter" does not refer to end-of-line (there is a separate property for that).
%     % Default textscan behaviour is to detect which of LF, CR, or CR+LF to use for end-of-line.
%     % delimiter='' probably means to ignore "
%     temp = textscan(fileId, '%s', 'delimiter', '', 'whitespace', '');
%     assert(numel(temp) == 1, 'textscan unexpectedly returned a non-scalar cell array.')
% 
%     rowsList = temp{1};
    
    fc = fread(fileId);
    rowsList = strsplit(char(fc'), linebreakRegexp, 'DelimiterType', 'RegularExpression', 'CollapseDelimiters', false)';
    


    fclose(fileId);
end
