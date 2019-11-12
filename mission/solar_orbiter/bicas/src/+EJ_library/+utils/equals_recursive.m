%
% Generic function for comparing MATLAB variables & data structures recursively.
%
% Can display whether variables are equal or not. If they differ, it will tell what the first found diffMsg is. Can
% handle numeric, char, logical, arrays, cells (incl. cell arrays), structs, objects not in arrays.
%
%
%
% ARGUMENTS
% =========
% a, b      : The data structures/variables that are to be compared.
% varargin  : Arguments as interpreted by EJ_library.utils.interpret_settings_args.
%   Possible settings 
%    .rowColVecEqual : True iff row and column vectors (1D arrays) should be considered equal despite having length in different dimensions/indices.
%    .nanEquals      : True iff NaN equals NaN.
%    .epsilon        : Permitted difference between numeric values.
%    Default values if field are not specified.
%
%
% RETURN VALUES
% =============
% result    : Whether the structs are equal or not.
%       0 : a and b differ
%       1 : a and b are equal
% diffLoc   : Location of difference between a and b.
% diffMsg   : Human readable text string describing the first difference found, if any. If no difference, then empty
%             string.
%
% Syntax for specifying location ("branching locations") in MATLAB variable is approximately
% the same as the MATLAB syntax:
%   -Component in array (including cell arrays with numel~=1): (i)
%   -Single-component cell array: {1}
%   -Field in struct: .<field name>
%
%
% USE CASES
% =========
% -Comparing similar cdfs (as represented as data structures; already loaded with other code), e.g. Solar Orbiter master CDFs.
% -Test code, comparing structs.
%
%
% NOTES
% =====
% NOTE: Empty strings 1x0 and 0x0 ('') count as unequal.
% NOTE: Very similar to isequal, isequaln, isequalwithequalnans (all recursive) except that this function
%    -Can pinpoint where one dissimilarity is.
%    -Has extra equality setting(s)
% NOTE: Distinguishes between different numeric types (counts as different), e.g. int64, double.
%
%
% Initially created 2016-07-0x by Erik P G Johansson.
%
function [result, diffLoc, diffMsg] = equals_recursive(a, b, varargin)
%
% QUESTION: Name?
%   PROPOSAL: equals
%   PROPOSAL: differ
%   PROPOSAL: compare
%
% PROPOSAL: Policies for different comparisons.
%     Display where diffMsg is.
%     Handles/objects : (1) Same handle (object), or (2) same classes and equal properties/fields, or (3) same properties/fields, count as equals.
%     The order of struct fields.
%     How to handle different classes (types) for numerical values: int64, in32, double and so on.
%     Approximate numeric equality, e.g. abs(a-b) < eps.
%
% PROPOSAL: Return list of all found differences, suiteable for reading on screen.
%    PROPOSAL: Return lists of all possible diffLoc & diffMsg.
% PROPOSAL: Implement exceptions for specific struct/object field names.
%    PROPOSAL: Specify by name.
%       PROPOSAL: settings.field_name_exceptions = {'fn1', 'fn2', ...}
%
% PROPOSAL: Change implementation to something that reuses code for all "branching".
%
% PROPOSAL: Function pointer to function for comparing numeric arrays.
%   PROPOSAL: Can customize comparison.
%       PRO: NaN, approximate equality
%
% PROPOSAL: Option for special treatment of string comparison. 0x0 and 1x0 could optionally count as equal.
%   NOTE: rowColVecEqual = true seems to have this effect for some reason. Bug?!
%
% TODO-NEED-INFO: Can objects be in array?!
% NOTE: There is something weird about objects.
% numel(containers.Map) == 1
% size(containers.Map) == [0, 1]
%
% numel(dataobj(...)) == 1
% size(dataobj(...)) == [1, 1]
%
% SPECIAL CASE: Empty arrays of different types. ==> Will not check for types if recursing over elements.



%==================================
% ASSERTIONS: Check argument types
%==================================
if ~numel(a) && ~isnumeric(a) && ~iscell(a) && ~ischar(a)   &&   ~numel(b) && ~isnumeric(b) && ~iscell(b) && ~ischar(b)
    error('Can not handle the type of variable(s)');
end



%======================
% Set default settings
%======================
DEFAULT_SETTINGS.rowColVecEqual = 0;
DEFAULT_SETTINGS.nanEquals      = 1;
DEFAULT_SETTINGS.epsilon        = 0.0;
%Settings = EJ_library.utils.add_struct_to_struct(settings, DEFAULT_SETTINGS, ...
%    struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'bothAreStructs', 'Error'));
Settings = EJ_library.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
EJ_library.utils.assert.struct2(Settings, {'rowColVecEqual', 'nanEquals', 'epsilon'}, {})   % ASSERTION



% Default return values.
result = 0;
diffLoc = '';
diffMsg = [];



if ~strcmp(class(a), class(b))
    % CASE: Class differs
    
    % NOTE: Must check equality of types first, since empty arrays will not be iterated over, and the comparison will
    % not be made checking the individual element.
    
    diffMsg = 'MATLAB classes (variable "types") differ.';
    result = false;
    return

elseif (numel(a) ~= 1) || (numel(b) ~= 1)
    % CASE: Same class, non-single-element arrays.
    
    a = canonicalize_size(a, Settings);
    b = canonicalize_size(b, Settings);

    if length(size(a)) ~= length(size(b))
        diffMsg = 'Number of dimensions differ.';
    elseif ~isequal(size(a), size(b))
        diffMsg = 'Array sizes differ.';
    else
        % CASE: a and b have identical array sizes.
        
        % This check is not necessary, but should speed up.
        if (isnumeric(a) && isnumeric(b)) || (ischar(a) && ischar(b))
             if ~equals_samesize_numeric_char_arrays(a, b, Settings)
                 diffMsg = 'Numeric/char arrays differ.'
             else
                 result = 1;
             end
        else
            for i = 1:numel(a)
                % NOTE: Technically set the values of diffLoc, diffMsg to the default values again.
                [result, diffLoc, diffMsg] = EJ_library.utils.equals_recursive(a(i), b(i), Settings);
                
                if result == 0
                    diffLoc = sprintf('(%i)%s', i, diffLoc);
                    return
                end
            end
            result = 1;
            %diffLoc = [];   % Must be reset?
        end
    end

elseif isobject(a) && isobject(b)
    % CASE: Neither a nor b is an array.
    
    if ~strcmp(class(a), class(b))     % Unnecessary?!!
        diffMsg = 'Class names differ.';
        return
    end
        
    fnList = fieldnames(a);
    for i = 1:length(fnList)
        fn = fnList{i};
        
        % NOTE: Technically set the values of diffLoc, diffMsg to the default values again.
        [result, diffLoc, diffMsg] = EJ_library.utils.equals_recursive(a.(fn), b.(fn));
        if result == 0
            diffLoc = sprintf('.%s%s', fn, diffLoc);
            return
        end
        
        % Comparing with MATLAB function, just to be sure the code can detect all differences.
        % This is useful for e.g. containers.Map.
        % MATLAB R2016a: "isequalwithequalnans is not recommended. Use ISEQUALN instead."
        result = isequalwithequalnans(a, b);
    end
    
elseif (isnumeric(a) && isnumeric(b)) || (ischar(a) && ischar(b)) || (islogical(a) && islogical(b))
    
    if ~equals_samesize_numeric_char_arrays(a, b, Settings)
        diffMsg = 'Numeric/char values differ.';
    else
        result = 1;
    end
        
elseif iscell(a) && iscell(b)
    
    [result, diffLoc, diffMsg] = EJ_library.utils.equals_recursive(a{1}, b{1}, Settings);   % RECURSIVE CALL
    if result == 0
        diffLoc = sprintf('{1}%s', diffLoc);
    end
    
elseif isstruct(a) && isstruct(b)
    
    % NOTE: length(setxor(fieldnames(a), fieldnames(b)))
    % does not test for the order of field names.
    mismatchingFieldNames = setxor(fieldnames(a), fieldnames(b));
    if ~isempty(mismatchingFieldNames)
        % IMPLEMENTATION NOTE: Important to write out fieldnames. Makes it easier to correct e.g. missing fieldnames.
        diffMsg = sprintf('Structs have different sets of field names. mismatchingFieldNames: %s', strjoin(mismatchingFieldNames, ', '));
    else
        fnList = fieldnames(a);
        for i = 1:length(fnList)
            fn = fnList{i};
            [result, diffLoc, diffMsg] = EJ_library.utils.equals_recursive(a.(fn), b.(fn), Settings);   % RECURSIVE CALL
            if result == 0
                diffLoc = sprintf('.%s%s', fn, diffLoc);
                return
            end
        end
    end

else
    diffMsg = 'Different types of variables, or types of variables the function can not handle.';    
end

end

%###################################################################################################

% Function for comparing same-size numeric/char arrays.
%
% Exists to provide one, single place to set the settings for what such a comparison means, e.g. NaN.
function result = equals_samesize_numeric_char_arrays(a, b, Settings)
    % NOTE: isequal interpret chars as numbers which can be equal to actual numbers, e.g.
    % isequal(65, 'A') == 1.
    % NOTE: isequal([], []) == 1
    % NOTE: isequal(NaN, NaN) == 0.
    % NOTE: isequaln(NaN, NaN) == 1. isequaln should otherwise be identical to isequal.
    % NOTE: isequalwithequalnans(NaN, NaN) == 1
    % NOTE: isequaln is not available in MATLAB 2009a.
    % NOTE: isequal, isqualn are recursive over structs, arrays and cell arrays.

    if isnumeric(a)
        if Settings.nanEquals
            %result = isequalwithequalnans(a, b);
            result = EJ_library.utils.approx_equals(a, b, Settings.epsilon, 'NaN equal to itself');
        else
            %result = isequal(a, b);
            result = EJ_library.utils.approx_equals(a, b, Settings.epsilon, 'NaN unequal to itself');
        end
    elseif ischar(a) || islogical(a)
        result = isequal(a, b);
    else
        error('a is neither numeric nor char.')
    end
    
end



function x = canonicalize_size(x, Settings)
    if Settings.rowColVecEqual && (numel(x) == length(x))
        % CASE: x is a 1D vector along an unknown dimension.
        x = reshape(x, [numel(x), 1]);   % Make into column vector.
    end    
end
