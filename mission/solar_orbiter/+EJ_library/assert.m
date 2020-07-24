%
% Class for static methods for creating assertions.
%
% NOTE: MATLAB itself (at least MATLAB R2009a) has a function "assert" which is useful for simpler assertions.
%
%
% POLICY
% ======
% Functions should be named as propositions (when including the class name "assert") which are true if the assertion
% function does not yield error. "castring" refers to arrays of char, not the concept of "strings" which begins with
% MATLAB 2017a and later.
%
%
% RATIONALE, REASONS FOR USING INSTEAD OF MATLAB's "assert"
% =========================================================
% PRO: Still more convenient and CLEARER for some cases.
%       PRO: More convenient for natural combinations of checks.
%           Ex: Structs (is struct + set of fields)
%           Ex: Strings (is char + size 1xN)
%           Ex: "String set" (cell array + only strings + unique strings)
%           Ex: "String list" (cell array + only strings)
%           Ex: Function handle (type + nargin + nargout)
%       PRO: Can have more tailored default error messages: Different messages for different checks within the same
%           assertions.
%           Ex: (1) Expected struct is not. (2) Struct has wrong set of fields + fields which are superfluous/missing.
% CON: Longer statement because of packages.
% PRO: MATLAB's assert requires argument to logical, not numerical (useful if using 0/1=false/true).
% --
% Ex: assert(strcmp(datasetType, 'DERIV1'))
%     vs EJ_library.assertions.castring_in_set(datasetType, {'DERIV1'})
% Ex: assert(any(strcmp(s, {'DERIV1', 'EDDER'})))
%     vs EJ_library.assertions.castring_in_set(datasetType, {'EDDER', 'DERIV1'})
% Ex: assert(isstruct(s) && isempty(setxor(fieldnames(s), {'a', 'b', 'c'})))
%     vs EJ_library.assertions.is_struct_w_fields(s, {'a', 'b', 'c'})
%
%
% NAMING CONVENTIONS
% ==================
% castring : Character Array (CA) string. String consisting 1xN (or 0x0?) matrix of char.
%            Name chosen to distinguish castrings from the MATLAB "string arrays" which were introduced in MATLAB R2017a.
%
%
% Initially created 2018-07-11 by Erik P G Johansson.
%
classdef assert
%
% TODO-DECISION: Use assertions on (assertion function) arguments internally?
% PROPOSAL: Add argument for name of argument so that can print better error messages.
% PROPOSAL: Optional error message (string) as last argument to ~every method.
%   CON: Can conflict with other string arguments.
%       Ex: Method "struct".
% PROPOSAL: Optional error message identifier as second-last argument to ~every method (see "error").
%   CON: Can conflict with other string arguments.
%
% PROPOSAL: Assertion for checking that multiple variables have same size in specific dimensions/indices.
%   See BOGIQ for method "size".
%
% PROPOSAL: Create class with collection of standardized non-trivial "condition functions", used by this "assert" class.
%           Use an analogous naming scheme.
%   PRO: Can use assertion methods for raising customized exceptions (not just assertion exceptions), e.g. for UI
%        errors.
%   PRO: Useful for creating more compact custom-made assertions.
%       Ex: assert(isscalar(padValue) || is_castring(padValue))
%       Ex: assert(<castring> || <struct>)
%       Ex: Checking settings values.
%           Ex: assert(ischar(defaultValue) || isnumeric(defaultValue) || is_castring_set(defaultValue))
%   PRO: Can use assertion methods for checking state/conditions (without raising errors).
%       Ex: if <castring> elseif <castring_set> else ... end
%   CON: Not clear what the conditions should be, and if they should some assertions themselves. Input checks or assume
%       the nominal condition
%       Ex: castring_set: Should the return value be whether the input argument is
%           A cell array of unique strings.
%           A cell array of unique strings (assertion: cell array)
%           A cell array of unique strings (assertion: cell array of strings) ???
%   PROPOSAL: Name "cond".
%   Ex: vector, struct
%   Ex: Because of comments?: dir_exists, file_exists
%   Ex: castring?
%   --
%   PROPOSAL: Have methods return a true/false value for assertion result. If a value is returned, then raise no assertion error.
%       PRO: Can be used with regular assert statements.
%           PRO: MATLAB's assert can be used for raising exception with customized error message.
%       CON: "assert" is a bad name for such a class.
%
% PROPOSAL: Static variable for error message identifier.
%   PRO: Can set differently in BICAS.
%
% PROPOSAL: Function for asserting that some arbitrary property is identical for an arbitrary set of variables.
%   Function handle defines the property: argument=variable, return value=value that should be identical for all
%   variables.
%   Ex: Size for some set of indices.
%       Ex: Range of first index (CDF Zvar records).
%           Ex: Can treat cell arrays specially: Check the components of cell array instead.
%
% PROPOSAL: Functions for asserting line breaks.
%   TODO-DECISION: Which set of functions.
%       PROPOSAL: Assert ending LF (assert not ending CR+LF).
%       PROPOSAL: Assert all linebreaks are LF (no CR+LF).
%       PROPOSAL: Assert all linebreaks are LF (no CR+LF). Require ending linebreak.
%
% PROPOSAL: Assert string sets equal
%   Ex: write_dataobj
%
% PROPOSAL: Utility function for error() and assert() that always use EJ_library.assert.ASSERTION_EMID.
%   CON: Not that much shorter.
%           error(EJ_library.assert.ASSERTION_EMID, msg)
%       ==> EJ_library.assert.error(msg)
%   CON: Gets method EJ_library.assert.assert() which is a bad name.
%       CON: Private function so matters less.
%   CON: Slower?
%       CON: Only applies to an internal assert function. Internal "error" function only triggered when assertion
%            alraedy has failed.
%   PRO: Easier to be consistent.
%   PRO: Clearer.
%       PRO: Fewer arguments.
%       



    properties(Constant)
        % EMID = Error Message ID
        %ERROR_EMID     = 'assert:IllformedAssertion'
        ASSERTION_EMID = 'assert:Assertion'
    end

    
    
    methods(Static)
        
        % NOTE: Empty string literal '' is 0x0.
        % NOTE: Accepts all empty char arrays, e.g. size 2x0 (2 rows).
        function castring(s)
            % PROPOSAL: Only accept empty char arrays of size 1x0 or 0x0.
            if ~ischar(s)
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected castring (0x0, 1xN char array) is not char.')
            elseif ~(isempty(s) || size(s, 1) == 1)
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected castring (0x0, 1xN char array) has wrong dimensions.')
            end
        end
        
        
        
        % Assert that ENTIRE string matches one of potentially many regexps.
        %
        % NOTE: Will permit empty strings to match a regular expression.
        %
        %
        % ARGUMENTS
        % =========
        % s      : String
        % regexp : (1) String. Regular expressions.
        %          (2) Cell array of strings. List of regular expressions.
        %              NOTE: Must be non-empty array.
        function castring_regexp(s, regexp)
            if ~any(EJ_library.str.regexpf(s, regexp))
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'String "%s" (in its entirety) does not match any of the specified regular expressions.', s)
            end
        end
        
        
        
        % Argument is a cell matrix of UNIQUE strings.
        function castring_set(s)
            % NOTE: Misleading name, since does not check for strings.
            
            if ~iscell(s)
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected cell array of unique strings, but is not cell array.')
                
            % IMPLEMENTATION NOTE: For cell arrays, "unique" requires the components to be strings. Therefore does not
            % check (again), since it is probably slow.
            elseif numel(unique(s)) ~= numel(s)
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected cell array of unique strings, but not all strings are unique.')
            end
        end
        
        
        
        function castring_sets_equal(set1, set2)
            % NOTE/BUG: Does not require sets to have internally unique strings.
            
            if ~isempty(setxor(set1, set2))
                error(EJ_library.assert.ASSERTION_EMID, 'The two string sets are not equivalent.')
            end
        end
        
        
        
        function castring_in_set(s, strSet)
        % PROPOSAL: Abolish
        %   PRO: Unnecessary since can use assert(ismember(s, strSet)).
        %       CON: This gives better error messages for string not being string, for string set not being string set.
        
            EJ_library.assert.castring_set(strSet)
            EJ_library.assert.castring(s)
            
            if ~ismember(s, strSet)
                error(EJ_library.assert.ASSERTION_EMID, 'Expected string in string set is not in set.')
            end
        end
        
        
        
        % NOTE: Can also be used for checking supersets.
        % NOTE: Both string sets and numeric sets
        function subset(strSubset, strSuperset)
            
            if ~EJ_library.utils.subset(strSubset, strSuperset)
                error(EJ_library.assert.ASSERTION_EMID, 'Expected subset/superset is not a subset/superset.')
            end
        end
        
        
        
        function scalar(x)
            if ~isscalar(x)
                error(EJ_library.assert.ASSERTION_EMID, 'Variable is not scalar as expected.')
            end
        end
        
        
        
        % Either regular file or symlink to regular file (i.e. not directory or symlink to directory).
        % Cf. path_is_available
        function file_exists(filePath)
            if ~(exist(filePath, 'file') == 2)
                error(EJ_library.assert.ASSERTION_EMID, 'Expected existing regular file (or symlink to regular file) "%s" can not be found.', filePath)
            end
        end
        
        
        
        function dir_exists(dirPath)
            if ~exist(dirPath, 'dir')
                error(EJ_library.assert.ASSERTION_EMID, 'Expected existing directory "%s" can not be found.', dirPath)
            end
        end
        
        
        
        % Assert that a path to a file/directory does not exist.
        %
        % Useful if one intends to write to a file (without overwriting).
        % Dose not assume that parent directory exists.
        function path_is_available(path)
        % PROPOSAL: Different name
        %   Ex: path_is_available
        %   Ex: file_dir_does_not_exist
        
            if exist(path, 'file')
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Path "%s" which was expected to point to nothing, actually points to a file/directory.', path)
            end
        end



        % ARGUMENTS
        % ==========
        % requiredFnSet : Cell array of required field names.
        % optionalFnSet : (a) Cell array of optional field names (i.e. allowed, but not required)
        %                 (b) String constant 'all' : All fieldnames are allowed but not required. This is likely only
        %                     meaningful when requiredFnSet is non-empty (not a requirement).
        %
        function struct(S, requiredFnSet, optionalFnSet)
            % PROPOSAL: Have it apply to a set of strings (e.g. fieldnames), not a struct as such.
            % PROPOSAL: Let optionalFnSet be optional (empty by default).
            %   PRO: Shorter for the most common case.
            %   PRO: Backward-compatibility with some of the syntax for struct(predecessor assertion function).
            %   CON: Bad for future extensions of function.
            %
            % PROPOSAL: Be able to (optionally) specify properties of individual fields.
            %   PROPOSAL: Arguments (one or many in a row) with prefix describe properties of previous field.
            %       Ex: 'fieldName', '-cell'
            %       Ex: 'fieldName', '-scalar'
            %       Ex: 'fieldName', '-double'
            %       Ex: 'fieldName', '-castring'
            %       Ex: 'fieldName', '-vector'
            %       Ex: 'fieldName', '-column vector'
            %       PRO: Can be combined with recursive scheme for structs, which can be regarded as an extension of
            %            this scheme. In that case, a cell array is implicitly interpreted as the assertion that the
            %            field is a struct with the specified (required and optional) subfields.
            %
            % PROPOSAL: Recursive structs field names.
            %   TODO-DECISION: How specify fieldnames? Can not use cell arrays recursively.
            %   PROPOSAL: Define other, separate assertion method.
            %   PROPOSAL: Tolerate/ignore that structs are array structs.
            %   PROPOSAL: struct(S, {'PointA.x', 'PointA.y'}, {'PointA.z'})
            %   PROPOSAL: struct(S, {'PointA', {'x', 'y'}}, {'PointA', {'z'}})   % Recursively
            %       Cell array means previous argument was the parent struct.
            %   PROPOSAL: struct(S, {'name', 'ReqPointA', {{'reqX', 'reqY'}, {'optZ'}}}, {'OptPointB', {{'reqX', 'reqY'}, {'optZ'}}})
            %       Cell array means previous argument was the parent struct.
            %       Groups together required and optional with every parent struct.
            %       PRO: Optional fields can be structs with both required and optional fields, recursively.
            %       PRO: Can be implemented recursively(?).
            %   TODO-NEED-INFO: Required & optional is well-defined?
            %   CON: Rarely needed.
            %       CON: Maybe not
            %           Ex: Settings structs
            %           Ex: EJ_library.so.abp.find_DSGs_by_filename_fields
            %   CON-PROPOSAL: Can manually call EJ_library.assert.struct multiple times, once for each substruct,
            %                 instead (if only required field names).
            %
            % PROPOSAL: Assertion: Intersection requiredFnSet-optionalFnSet is empty.
            
            structFnSet          = fieldnames(S);
            
            missingRequiredFnSet = setdiff(requiredFnSet, structFnSet);
            
            % disallowedFnSet = ...
            if iscell(optionalFnSet)
                disallowedFnSet = setdiff(setdiff(structFnSet, requiredFnSet), optionalFnSet);
            elseif isequal(optionalFnSet, 'all')
                disallowedFnSet = {};
            else
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Illegal optionalFnSet argument. Is neither cell array or string constant "all".')
            end
            
            % Give error, with an actually useful error message.
            if ~isempty(missingRequiredFnSet) || ~isempty(disallowedFnSet)
                missingRequiredFnListStr = strjoin(missingRequiredFnSet, ', ');
                disallowedFnListStr      = strjoin(disallowedFnSet,      ', ');

                error(EJ_library.assert.ASSERTION_EMID, ...
                    ['Expected struct has the wrong set of fields.', ...
                    '\n    Missing fields:           %s', ...
                    '\n    Extra (forbidden) fields: %s'], missingRequiredFnListStr, disallowedFnListStr)
            end
        end



        % ARGUMENTS
        % =========
        % nArgin   : Value returned by nargin() for function handle.
        %            abs(nArgin)  = Number of arguments, where varargin counts as one.
        %            sign(nArgin) = 0 or 1: No varargin.
        %                          -1     : Has varargin.
        % nArgout : Value returned by nargout() for function handle.
        %            Analogous to nArgin.
        % 
        % NOTE: Can not handle function handle that does not point to existing function(!)
        %
        % NOTE: Method's USEFULNESS IS LIMITED due to below reasons:
        %   NOTE: Can not distinguish between functions that 
        %       (1) only accept exactly N arguments, and
        %       (2) accept a variable-number of arguments, including N arguments (using varargin).
        %       Analogous problem with nargout.
        %   NOTE: Empirically, nargout(anonymousFunction) == -1, always.
        %
        function func(funcHandle, nArgin, nArgout)
            % PROPOSAL: Make nArgin/nArgout able to simultaneously accept multiple nargin/nargout values.
            %   Ex: Accept nargin  = 3, -1,-2,-3 (accept all functions that seemingly can accept three arguments)
            %   Ex: Accept nargout = 3, -1,-2,-3 (accept all functions that seemingly can return three values)
            %   PROPOSAL: Submit exact set (numeric array) of accepted values.
            
            if ~isa(funcHandle, 'function_handle')
                error(EJ_library.assert.ASSERTION_EMID, 'Expected function handle is not a function handle.')
            end
            if nargin(funcHandle) ~= nArgin
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected function handle ("%s") has the wrong number of input arguments. nargin()=%i, nArgin=%i', ...
                    func2str(funcHandle), nargin(funcHandle), nArgin)
            elseif nargout(funcHandle) ~= nArgout
                % NOTE: MATLAB actually uses term "output arguments".
                error(EJ_library.assert.ASSERTION_EMID, ...
                    ['Expected function handle ("%s") has the wrong number of output', ...
                    ' arguments (return values). nargout()=%i, nArgout=%i'], ...
                    func2str(funcHandle), nargout(funcHandle), nArgout)
            end
        end



        function isa(v, className)
            if ~isa(v, className)
                error(EJ_library.assert.ASSERTION_EMID, 'Expected class=%s but found class=%s.', className, class(v))
            end
        end
        
        

        % Assert v has a non-one size in at most one dimension.
        % NOTE: Excludes 0x0, e.g. [] and ''.
        %
        % NOTE: MATLAB's "isvector" function uses different criterion which excludes length in third or higher
        % dimension.
        %   isvector(ones(0,1)) == 1
        %   isvector(ones(1,0)) == 1
        %   isvector(ones(1,1,3)) == 0
        %   isvector(ones(1,1,0)) == 0
        %
        function vector(v)
            % PROPOSAL: Optional extra argument that specifies the length.
            % PROPOSAL: Better name.
            %   PRO: "vector" is ambiguous w.r.t. number of dimensions.
            %   PROPOSAL: Name that implies 1D vector: vector_1D, vec_1D
            %   PROPOSAL: Name that implies exactly one non-zero size dimension.
            %       PROPOSAL: true_vector, strict_vector, strict_vec_1D, vec_strict_1D
            % PROPOSAL: Permit 0x0.
            %   PROPOSAL: Separate constant to permit. '0x0', 'permit 0x0'.
            %   PROPOSAL: Separate method to permit.
            %       PROPOSAL: common_vec_1D, nonstrict_vec_1D.
            
%             dims = size(v);
%             dims(dims==1) = [];
%             if numel(dims) > 1
            if sum(size(v) ~= 1) > 1
                sizeStr = sprintf('%ix', size(v));
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected vector, but found variable of size %s.', sizeStr(1:end-1))
            end
        end
        
        
        
        % Assert that v has a specific size, in all or some dimensions/indices.
        % NOTE: Should eventually be replaced by sizes()?
        %
        % ARGUMENTS
        % =========
        % v               : Variable which size will be asserted.
        % sizeConstraints : 1D vector with the sizes of the corresponding indices/dimensions. A component value of NaN
        %                   means that the size of that particular dimension will not be checked. Higher dimensions
        %                   which are not specified are implicitly one.
        %
        function size(v, sizeConstraints)
            % PROPOSAL: Apply the same size constraint to an arbitrary number of variables.
            %
            % Cf EJ_library.assert.sizes.
            
            
            % ASSERTION
            EJ_library.assert.vector(sizeConstraints)
            
            sizeV = size(v);
            
            % Enforce column vectors.
            sizeV           = sizeV(:);
            sizeConstraints = sizeConstraints(:);
            
            nSizeV           = numel(sizeV);
            nSizeConstraints = numel(sizeConstraints);
            
            % Enforce that sizeV and sizeConstraints have equal size by adding components equal to one (1).
            % NOTE: MATLAB's "size" function always returns at least a 1x2 vector.
            if (nSizeV < nSizeConstraints)
                sizeV           = [sizeV;           ones(nSizeConstraints-nSizeV, 1)];
            else
                sizeConstraints = [sizeConstraints; ones(nSizeV-nSizeConstraints, 1)];
            end
            
            % Overwrite NaN values with the actual size values for those indices.
            bIgnore = isnan(sizeConstraints);
            sizeConstraints(bIgnore) = sizeV(bIgnore);

            % ASSERTION: The actual assertion
            assert( all(sizeV == sizeConstraints), ...
                EJ_library.assert.ASSERTION_EMID, 'Variable does not have the expected size.')
        end



        % Simultaneously
        % ** check the size of multiple variables,
        % ** whether specified arbitrary-sizes dimensions are identical across dimensions and variables
        % ** return above arbitarily-sized dimensions, so that the caller can e.g. check them.
        %
        % Ex: Check that variables representing zVariables have the right non-record dimensions, and the same number of
        % records (size in first dimension). Return number of records.
        %
        %
        % ARGUMENTS
        % =========
        % varargin: Arbitrary number of argument pairs below:
        %   var            : variable value which size will be tested.
        %   sizeConstraint : 1D vector with integers specifying size of corresponding argument "var".
        %       Negative integer : Arbitrary dimension size which must match between all arguments "var".
        %                          Must be numbered -1, -2, ... , -N
        %       NaN              : Arbitrary dimension size independent of other dimensions. No corresponding return
        %                          value.
        % 
        %
        % RETURN VALUES
        % =============
        % Size of dimensions labelled with negative integers, in order -1, -2, ... .
        %
        %
        % NOTES
        % =====
        % FINISHED (except assertions on input) but uncertain if
        %   has sensible design
        %   has sensible name
        %   should replace EJ_library.assert.size().
        % NOTE: "Backward-compatible" with EJ_library.assert.size().
        % NOTE: Returning values in principle makes it an assertion+functionality.
        % NOTE: One needs to add semicolon to end of row, since has return values.
        %
        function [varargout] = sizes(varargin)
            % PROPOSAL: Be able to separate size constraints to multiple variables, but specify that certain indices
            %           have to be identical in size (but arbitrary) between variables.
            %
            %   PROPOSAL: Somehow be able to state that a variable is a 1D vector, regardless of which index is not size one.
            %       PROPOSAL: sizeConstraints = {N}, one numeric value (N, negativeValue, NaN).
            %       PROPOSAL: sizeConstraints = {'1D vector', N}
            %       PROPOSAL: Prepend sizeConstraints argument with string constant "vector", "1D", "1D vector".
            %   ~CON/NOTE: Can not assert equal size for variables with arbitrary number of dimensions.
            %
            %   PROPOSAL: Same dimensions in all dimensions except those specified.
            %       PRO: Better handling of "high dimensions to infinity".
            %
            % PROPOSAL: Count all negative numbers as sizes, not just successive integers.
            %           Return values are negative integers in order.
            % PROPOSAL: Abolish NaN (force caller to use negative integer size).
            %   PRO: Simplifies implementation.
            %   PRO: Saves NaN for possible future special value for other functionality.
            %   CON: Clearer that NaN refers to arbitrary value.
            %
            % PROPOSAL: Use method to replace EJ_library.assert.size().
            %   PROPOSAL: Keep name "sizes".
            %       CON: Bad name when testing just one variable.
            %   PROPOSAL: Use (take over) name "size".
            %       PRO: Can see constraining many variables as an extension of the base functionality.
            
            nArgs = numel(varargin);

            sizeArray            = [];
            sizeConstrArray = [];

            %=====================================================================================
            % Read arguments and store values in a data structure that is suitable for algorithm.
            %=====================================================================================
            iArgPair = 0;    % NOTE: iArgPair=1 represents the first argument pair.
            while true
                %====================
                % Read argument pair
                %====================
                if nArgs == 2*iArgPair
                    break
                elseif nArgs >= 2*(iArgPair+1)
                    iArgPair = iArgPair + 1;
                else
                    error(EJ_library.assert.ASSERTION_EMID, 'Ill-formed assertion. Number of arguments is not even.')
                end

                sizeArg           = size(varargin{2*iArgPair-1});
                sizeConstraintArg =      varargin{2*iArgPair  };
                EJ_library.assert.vector(sizeConstraintArg)
                
                %====================================================
                % Add argument pair values into same-sized 1D arrays
                %====================================================                
                % Force column arrays.
                sizeArg           = sizeArg(:);
                sizeConstraintArg = sizeConstraintArg(:);                
                
                % Pad the smallest arrays with ones until both have same size
                % -----------------------------------------------------------
                % NOTE: padarray pads in first dimension by default.
                %       'post' : Pads after the last array element along each dimension.
                nDiff = numel(sizeArg) - numel(sizeConstraintArg);
                if nDiff >= 1   sizeConstraintArg = padarray(sizeConstraintArg,  nDiff, 1, 'post');
                else            sizeArg           = padarray(sizeArg,           -nDiff, 1, 'post');
                end
                
                sizeArray       = [sizeArray;            sizeArg];
                sizeConstrArray = [sizeConstrArray; sizeConstraintArg];
            end
            
            %==============================================
            % ASSERTION: Explicitly stated dimension sizes
            %==============================================
            b = (sizeConstrArray >= 0);    % NOTE: (NaN >= 0) returns false.
            assert(all(sizeConstrArray(b) == sizeArray(b)), ...
                EJ_library.assert.ASSERTION_EMID, ...
                'At least one variable dimension size differs from an explicitly specified dimension size.')
            sizeConstrArray(b) = NaN;    % Effectively remove already checked size constraints from later checks.

            %=====================================================================================
            % ASSERTION: Linked dimension sizes must match each other (negative size constraints)
            %=====================================================================================
            sizeConstrSpecialValue = -1;
            while true
                b = (sizeConstrArray == sizeConstrSpecialValue);
                if ~any(b)
                    % CASE: sizeConstrSpecialValue is not used as a special value.
                    break
                end

                uniqueSizes  = unique(sizeArray(b));
                nUniqueSizes = numel(uniqueSizes);

                assert(nUniqueSizes == 1, ...
                    EJ_library.assert.ASSERTION_EMID, ...
                    'One or several variables have different sizes for dimensions labelled %i.', sizeConstrSpecialValue)

                % Assign return values.
                varargout{-sizeConstrSpecialValue} = uniqueSizes;

                sizeConstrArray(b) = NaN;    % Effectively remove already checked size constraints from later checks.
                
                sizeConstrSpecialValue = sizeConstrSpecialValue - 1;
            end
            
            %===========================================
            % ASSERTION: Only size constraint NaN left.
            %===========================================
            assert(all(isnan(sizeConstrArray)), ...
                EJ_library.assert.ASSERTION_EMID, ...
                ['Size constraints contains negative numbers that can not (are not supposed to)', ...
                ' be interpreted as constraints.'])

            %==============================================================================
            % NOTE: Ignore size constraint NaN. Does not need to be explicitly checked.
            %==============================================================================
        end



        % Assert that all values in a matrix are identical. Useful for e.g. checking that sizes of vectors are
        % identical.
        %
        %
        % ARGUMENT
        % ========
        % v : One of below:
        %     - matrix of numbers
        %     - matrix of characters
        %     - cell array of strings
        %     Does not work on:
        %     - cell array of numbers
        %     NOTE: Empty matrices are accepted.
        %
        function all_equal(v)
            % TODO-DECISION: How handle NaN?
            %   PROPOSAL: Count NaN as equal to itself.
            
            nUniques = numel(unique(v(:)));    % NOTE: Make 1D vector.
            nTotal   = numel(v);
            
            if (nUniques ~= 1) && (nTotal >= 1)
                error(EJ_library.assert.ASSERTION_EMID, ...
                    'Expected vector of identical values, but found %i unique values out of a total of %i values.', ...
                    nUniques, nTotal)
            end
        end
        
    end    % methods
end    % classdef
