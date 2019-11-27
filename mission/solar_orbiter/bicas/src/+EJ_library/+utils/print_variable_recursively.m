%
% Print the contents of a MATLAB variable recursively. Intended for getting an overview of the contents of complex
% variables consisting of nested structures, nested cell arrays, etc. useful for debugging.
%
%
% ARGUMENTS
% =========
% varName  : String that will be used as a (top-level) variable name.
% v        : The variable which contents will be printed.
% varargin : Settings as interpreted by EJ_library.utils.interpret_settings. See implementation.
%
%
% DISPLAYING STRINGS
% ==================
% The function displays every string in one of two different ways:
% (1) BAE = Begin After Equals : More compact. Better for strings with few/no linebreaks.
% (2) BNR = Begin on New Row   : Useful for strings with many linebreaks.
%   Setting stringsBnrMinNonemptyRows determines when to use which.
%
%
% NOTE: NOT COMPLETE. More cases to be added as they are needed.
% NOTE: Option maxRecursionDepth implemented but only crudely.
% NOTE: Will only print public properties for objects.
%
%
% First created by Erik P G Johansson 2017-03-21.
%
function print_variable_recursively(varName, v, varargin)
% PROPOSAL: Returna multi-rad sträng istf att skriva ut.
%   PRO: Can indent further.
%   PRO: Kan söka igenom (efter fältnamn).
%   PRO: Kan implementera automatiskt testkod.
%   PROPOSAL: Returnera som sträng om anroparen anger en returvariabel.
%
% PROPOSAL: Rename to something shorter.
%   PROPOSAL: print_vr
%
% PROPOSAL: Settings.maxRecursionDepth
%   TODO-DECISION: Vad skriva ut när når gränsen och alltså valt att inte skriva ut delkomponenter?
%       Ex: Når rekursionsgränsen när når komponent i struct array, dvs. t.ex. "s.Groups(2).Datasets(2)"? Vill inte
%       upprepa en massa identiska utskrifter (varje struct array-komponent har samma lista fält). Inte befogat med
%       någon sorts sammanfattning som är identisk för varje komponent. För cell arrays kan komponenterna vara olika så
%       där kan någon sammanfattning vara befogad.
%
% PROPOSAL: Line up equal signs among children of same parent.
%   PRO: More easily readable.
%   PRO: More elegant. ==> Can use function for more "formal" printouts.
%       Ex: EJ_library.atest.automatically_test_function, incl. exceptions.
%
% BUG: maxRecursionDepth<Inf + prinParentsSeparately=false will not print children with children (non-leafs).
% BUG? Settings.printParentSeparately does not work?
%
%
% Overview of cases
% -----------------
% NOT: All variables have an array size, even though one may to treat 1x1 as a special case.
% Cells
% Char: row vector counts as string
% Struct, object (class object)
% Numeric
% Function handle
%
% PROPOSAL: Criterion for displaying strings over multiple lines?
%   PROPOSAL: Min number of non-empty rows
%   PROPOSAL: Min number of non-linebreak chars
%   PROPOSAL: Multiple criteria. Use MLD when at least one of them is satisfied.
%
% TODO-DECISION: How handle strings (1xN char arrays) and general char arrays (2-D, N-D)?
%   TODO-DECISION: How handle strings (1xN char arrays)
%                  How handle long strings ~without linebreaks?
%                  How handle long strings with many linebreaks?
%       PROPOSAL: Separate printing mode which divides string into separate substrings. Return cell array of substrings.
%                 Caller can print strings on separate substrings, with indentation, prefix, everyone separately quoted.
%           PROPOSAL: Separate modes/reasons/criteria for splitting into separate substring.
%               NOTE: If splitting for other reasons than ~linebreak, then complete string is ambiguous if not using escape
%                     codes.
%                   NOTE: String is ALWAYS ambiguous if not using escape codes, e.g. tab.
%               NOTE: If splitting string into substrings which should themselves not be linebroken when printed, then
%                     line-breaking characters must be removed manually (even if not dispaly escape codes).
%                   PROPOSAL: Print with special function which tries to ignore (remove) linebreaking characters.
%               --
%               PROPOSAL: Split on linebreak.
%               PROPOSAL: Max substring length.
%                   TODO-DECISION: Before/after escape codes? How handle after escape codes?
%               PROPOSAL: Print comment with reason for splitting.
%                   CON: Must be returned, and printed by caller.
%
% PROPOSAL: Print multirow valueDisplayStrings using indentation for non-first row.
% PROPOSAL: Make "=" line up, at least among siblings.
% PROPOSAL: Rename "printParentSeparately" --> "printParentsSeparately" (plural)

% Define default settings.
DEFAULT_SETTINGS = [];
DEFAULT_SETTINGS.maxNArrayComponents       = Inf;      % Max number of array components to print.
DEFAULT_SETTINGS.stringsMaxDisplayLen      = Inf;      % Max length of printed strings.
DEFAULT_SETTINGS.maxRecursionDepth         = Inf;      % First function call is level 0. Counting may not be perfectly implemented.
DEFAULT_SETTINGS.indent                    = false;    % Visualisera underträd mha indentering. Implicerar en extra rad för föräldern till varje underträd.
DEFAULT_SETTINGS.indentationLength         = 4;        % Indentation length that is added per recursion level, if indentation is used.
DEFAULT_SETTINGS.printParentSeparately     = true;
DEFAULT_SETTINGS.stringsEscape             = true;     % Whether to display (some) special characters using escape codes.
%DEFAULT_SETTINGS.stringsBnrMinNonemptyRows = 3';       % Criterion for whether to use BNR.
DEFAULT_SETTINGS.stringsSsMaxLength        = 120;      % Max length on (pre-escaped) substrings before dividing into more lines.

Settings = EJ_library.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
EJ_library.utils.assert.struct2(Settings, fieldnames(DEFAULT_SETTINGS), {})    % Require DEFAULT_SETTINGS to contain all possible settings.

print_NESTED('', varName, v, 0, Settings)

end



% recursionDepth is used for determining the indentation.
function print_NESTED(fullParentName, varName, v, recursionDepth, Settings)
    if Settings.maxRecursionDepth < recursionDepth
        return
    end
    
    [canHaveChildren, childrenVList, childrenNamesList, valueDisplayStr] = interpret(v, Settings);
    

    
    if ~canHaveChildren
        % CASE: Node can not have children (under any circumstance).
        if Settings.indent
            prefixStr = indentation_str(recursionDepth * Settings.indentationLength);
        else
            prefixStr = fullParentName;
        end
        
        %============
        % Print leaf
        %============
        % NOTE: Arbitrarily chosen string width to make = line up, normally.
        fprintf(1, '%s%-20s = %s\n', prefixStr, varName, valueDisplayStr);
        
    else
        % CASE: Node may have children (but does not have to have this time).
        childrensFullParentName = [fullParentName, varName];
        
        %======================
        % Print parent (maybe)
        %======================
        if Settings.indent
            prefixStr = indentation_str(recursionDepth * Settings.indentationLength);
        else
            prefixStr = fullParentName;
        end        
        %if Settings.indent || Settings.printParentSeparately   % Why checking for .indent?
        if isempty(childrenVList) || Settings.printParentSeparately
            % NOTE: Want to always write parent if it has no children. Will oterwise not be reprented at all.
            %   Ex: Empty struct.
            fprintf(1, '%s%s   (size %s, class "%s")\n', prefixStr, varName, size_str(v), class(v));
        end
        
        %============================
        % Print children (recursive)
        %============================
        nChildren = length(childrenVList);
        if isstruct(v) && numel(v) == 1     % NOTE: Should not trigger for struct arrays
            nPrintComponents = nChildren;
        else
            nPrintComponents = min(length(childrenVList), Settings.maxNArrayComponents);
        end
        for i = 1:nPrintComponents
            print_NESTED(childrensFullParentName, childrenNamesList{i}, childrenVList{i}, recursionDepth+1, Settings)    % RECURSIVE CALL
        end
        if nPrintComponents < nChildren
            fprintf(1, '%s ... (too many array components; total %i components/cells)\n', ...
                [fullParentName, varName, childrenNamesList{i+1}], ...
                nChildren ...
                );
        end
    end
end



% Take a variable and either
% (1) create the appropriate one-row display string, or
% (2) identify and return the variable's children and their "relative names" (cf "relative paths").
%
% This function is NOT recursive, but is intended to be used by other recursive code.
% This means that this function decides which variable values should be interpreted as having children or not, e.g.
% character strings (which are really row vectors).
%
% RETURN VALUES
% =============
% Returns either
% (1)
%   canHaveChildren   = false
%   childrenVList     = {}
%   childrenNamesList = {}
%   displayValue      = '...'    % String, mostly one row.
% or (2)
%   canHaveChildren   = true
%   childrenVList     = {...}
%   childrenNamesList = {...}
%   displayValue      = []
%
%
% NOTE: DOES NOT IMPLEMENT ALL CASES.
% NOTE: canHaveChildren is useful for some types of printing when the "parent node" should be printed.
% NOTE: canHaveChildren can be true at the same time as there are no children.
%   Ex: Struct without fields, empty array/cell array.
% NOTE: Not entirely clear how to treat all cases
%   Ex: Empty char/numeric/cell array, empty struct -- Should it be printed on its own line even when there is no header line?
%   Ex: Single component char/numeric array         -- Should it be printed with an index   even when there is no header line?
%
% ~BUG: Single-component (non-array) structs will be displayed as a length-one struct array.
% BUG: Can not handle non-row, non-empty char arrays.

function [canHaveChildren, childrenVList, childrenNamesList, valueDisplayStr] = interpret(v, Settings)
% PROPOSAL: Replace "canHaveChildren" and "displayValue" with headerStr (always returned). headerStr can be used (1) as a
% summary for variables which can have children but do not, (2) for header/title lines when indenting, (3) optional .
%
    LF = char(10);
    
    childrenVList     = {};
    childrenNamesList = {};
    valueDisplayStr   = [];
    
    % IMPLEMENTATION NOTE: Cases are sorted in a hierarchical manner. The top-most category is number of
    % components/children (char arrays since row arrays and empty arrays need to be treated specially).
    % Second-top-most category is type of data.
    
    if ischar(v) && (is_row(v) || all(size(v) == [0,0]))
        % CASE: Regular string (one row or empty)
        
        if isempty(v)
            %====================
            % CASE: Empty string
            %====================
            canHaveChildren = false;    % A bit wrong.

            % Print empty string and array size. 
            % NOTE: Size 0x0 for regular empty strings, '', otherwise 1x0.
            valueDisplayStr = sprintf('''''   (size %s, class "%s")', size_str(v), class(v));
        else
            %========================
            % CASE: Non-empty string
            %========================
            canHaveChildren = false;
            
            % Optionally shorten content string
            if length(v) > Settings.stringsMaxDisplayLen
                valueContentStr = sprintf('%s...', v(1:Settings.stringsMaxDisplayLen));   % NOTE: Technically adding characters.
                valueCommentStr = sprintf('too many chars; first %i chars of %i', ...
                    Settings.stringsMaxDisplayLen, length(v));
            else
                valueContentStr = v;
                valueCommentStr = sprintf('%i chars', numel(valueContentStr));
            end
            
            % Split string into multiple substrings, to display on separate rows.
            valueDisplayStrList = split_content_str(valueContentStr, Settings.stringsSsMaxLength);
            
            if Settings.stringsEscape
                valueDisplayStrList = special_chars_2_escape_codes(valueDisplayStrList);
            end
            
            % Construct one single long string for display.
            valueDisplayStrList      = remove_linebreaks(valueDisplayStrList);                                               % Remove content linebreaks (not escape codes).
            valueDisplayStrList      = cellfun(@(x) (sprintf('''%s''', x)), valueDisplayStrList, 'UniformOutput', false);    % Quote every substring.
            valueDisplayStrList{end} = sprintf('%s   (%s)', valueDisplayStrList{end}, valueCommentStr);                      % Add comment string to last string.
            
            if numel(valueDisplayStrList) == 1
                valueDisplayStr = valueDisplayStrList{1};
            else
                % Have actual string substrings be displayed on separate rows.
                valueDisplayStr = [LF, strjoin(valueDisplayStrList, LF)];
            end

        end

    elseif numel(v) == 0
        % CASE: Empty matrix (char, non-char)
        
        if isnumeric(v) || islogical(v)
            canHaveChildren = false;    % A bit wrong.
            valueDisplayStr = sprintf('[]   (size %s, class "%s")', size_str(v), class(v));   % Print array size.
            
        elseif iscell(v)
            canHaveChildren = false;    % A bit wrong.
            valueDisplayStr = sprintf('{}   (size %s, class "%s")', size_str(v), class(v));   % Print array size.            
            
        elseif isstruct(v) || isobject(v)
            % NOTE: Only gets the PUBLIC properties of objects/classes.
% 
            canHaveChildren = false;
            valueDisplayStr = sprintf('(size %s, class "%s", with fields: %s)', size_str(v), class(v), strjoin(fieldnames(v), ', '));

        else
            error('Can not handle this variable size or type. class(v)="%s"', class(v))
        end
        
    elseif numel(v) == 1
        % CASE: 1x1

        if isnumeric(v)
            canHaveChildren = false;
            if imag(v) == 0
                % CASE: Real number
                valueDisplayStr = sprintf('%g', v);
            else
                % CASE: Complex number
                % NOTE: Want operator * in case imaginary part is NaN.
                valueDisplayStr = sprintf('%g + %g*i', real(v), imag(v));
            end
            
        elseif islogical(v)
            canHaveChildren = false;
            if v
                valueDisplayStr = 'true';
            else
                valueDisplayStr = 'false';
            end
 
        elseif isstruct(v) || isobject(v)
            % NOTE: Only gets the PUBLIC properties of objects/classes.

            canHaveChildren = true;
            fieldNamesList  = fieldnames(v);
            for iField = 1:length(fieldNamesList)
                childrenNamesList{iField} = ['.', fieldNamesList{iField}];
                childrenVList{iField}     = v.(fieldNamesList{iField});
            end

        elseif iscell(v)
            canHaveChildren      = true;
            childrenVList{1}     = v{1};
            childrenNamesList{1} = '{1}';

        elseif isa(v, 'function_handle')
            canHaveChildren = false;
            
            %valueDisplayStr = sprintf('%s: nargin=%i, nargout=%i (function handle)', func2str(v), nargin(v), nargout(v));
            valueDisplayStr = sprintf('%s (function handle)', func2str(v));

        else
            error('Can not handle this variable size or type.')
        end
        
    else
        % CASE: Non-empty, non-single element array
        
        if iscell(v)            
            canHaveChildren = true;
            for i = 1:numel(v)
                childrenVList{i}     = v{i};
                %childrenNamesList{i} = sprintf('{%i}', i);
                childrenNamesList{i} = sprintf('{%s}', index_str(size(v), i));
            end

        elseif isstruct(v) || isnumeric(v) || islogical(v) || isobject(v)
            canHaveChildren = true;
            for i = 1:numel(v)
                childrenVList{i}     = v(i);
                childrenNamesList{i} = sprintf('(%s)', index_str(size(v), i));
                %childrenNamesList{i} = sprintf('(%i)', i);
            end

        else            
            error('Can not handle this variable size or type (MATLAB class).')
        end
        
    end
end   % function



% Convert a "linear index" into "subscript values" (see "ind2sub"), i.e. converts a single linear matrix index (for matrix of arbitrary dimension) into
% multiple indices, one per matrix dimension.
function s = index_str(siz, i)
    iCellArray = cell(1, numel(siz), 1);
    [iCellArray{:}] = ind2sub(siz, i);
    s = regexprep(num2str([iCellArray{:}]), ' *', ',');
end



% Return display string representing the array size (the result of "size", i.e. all arrays incl. struct arrays, cell
% arrays) of an arbitrary variable, e.g. "3x4", "1x0".
function str = size_str(v)
    str = regexprep(num2str(size(v)), ' *', 'x');
end



% NOTE: MATLAB builtin "isrow" does not exist in MATLAB2009a.
function result = is_row(v)
    result = (ndims(v) == 2) && (size(v,1) == 1);
end



function s = indentation_str(N)
    s = repmat(' ', [1, N]);
end



% Remove characters that break lines when printed (with fprintf etc).
%
% ARGUMENTS
% =========
% s : String, or cell array of strings.
function s = remove_linebreaks(s)
    s = strrep(s, char(10), '');
    s = strrep(s, char(13), '');
end



% Replace some special characters with escape codes.
% Inefficient?
%
% ARGUMENTS
% =========
% s : String, or cell array of strings.
function s = special_chars_2_escape_codes(s)
    % PROPOSAL: Make into generic function.
    %   PRO: Automatic testing.
    
    % Escape codes for characters that will be displayed as such. Escape codes here are defined by sprintf.
    % \\ must come first. Will otherwise affect other substitutions.
    ESC_STR_LIST = {'\\', '\n', '\r', '\t'};

    for i = 1:numel(ESC_STR_LIST)
        escStr  = ESC_STR_LIST{i};
        charStr = sprintf(escStr);
        
        s = strrep(s, charStr, escStr);
    end
end



% Split string for display (if deemed necessary).
% Must be a string without escape codes.
function strList = split_content_str(s, maxSsLength)
    
    %LF = char(10);
    
    assert(maxSsLength>=1)
    
    % IMPLEMENTATION NOTE: Can not use regexp(s, '\n', 'split')
    % since we want to keep the linebreaking character in the substrings.
    % Split into substrings, where each substring contains the longest possible sequence of non-LF characters followed
    % by either (1) LF, or (2) end-of-string.
    %strList1 = regexp(s, '[^\n]*(\n|$)', 'split');
    
    % Split into substrings, where each substring contains the longest possible sequence of non-LF characters (that is
    % no longer than maxSsLength-1 characters) followed by exactly one character of either (1) a non-LF charcter (2) a
    % LF, or (3) end-of-string.
    strList = regexp(s, sprintf('[^\n]{0,%i}([^\n]|\n|$)', maxSsLength-1), 'match');

    % ASSERTION: Length of all substrings == length of original string
    assert(sum(cellfun(@length, strList)) == length(s))
end



% % Function which decides/determines whether string should be displayed using BNR.
% function useBnr = use_multiline_display(s, minNonemptyRows)
%     %EJ_library.utils.assert.castring(s)    % DEBUG
%     
%     strRowsList = regexp(s, '\n', 'split');
%     nNonemptyRows = numel(cellfun(@isempty, strRowsList));
%     useBnr = (nNonemptyRows >= minNonemptyRows);
% 
% end
