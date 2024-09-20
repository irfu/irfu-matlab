%
% Class for static methods for creating assertions.
%
% NOTE: MATLAB itself (at least MATLAB R2009a) has a function "assert" which is
% useful for simpler assertions.
%
%
% POLICY
% ======
% Functions should be named as propositions (when including the class name
% "assert") which are true if the assertion function does not yield error.
% "castring" refers to arrays of char, not the concept of "strings" which begins
% with MATLAB 2017a and later.
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
%       PRO: Can have more tailored default error messages: Different messages
%            for different checks within the same assertions.
%           Ex: (1) Expected struct is not. (2) Struct has wrong set of fields +
%               fields which are superfluous/missing.
% CON: Longer statement because of packages.
% PRO: MATLAB's assert requires argument to logical, not numerical (useful if
%      using 0/1=false/true).
% --
% Ex: assert(strcmp(datasetType, 'DERIV1'))
%     vs irf.assertions.castring_in_set(datasetType, {'DERIV1'})
% Ex: assert(any(strcmp(s, {'DERIV1', 'EDDER'})))
%     vs irf.assertions.castring_in_set(datasetType, {'EDDER', 'DERIV1'})
% Ex: assert(isstruct(s) && isempty(setxor(fieldnames(s), {'a', 'b', 'c'})))
%     vs irf.assertions.is_struct_w_fields(s, {'a', 'b', 'c'})
%
%
% NAMING CONVENTIONS
% ==================
% castring : Character Array (CA) string. String consisting 1xN (or 0x0?) matrix
%            of char. Name chosen to distinguish castrings from the MATLAB
%            "string arrays" which were introduced in MATLAB R2017a.
%
%
% Initially created 2018-07-11 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
classdef assert
  %
  % TODO-DEC: Use assertions on (assertion function) arguments internally?
  %   NOTE: Potentially slower.
  %
  % PROPOSAL: Internal utility function for error() and assert() that always use irf.assert.ASSERTION_EMID.
  %   CON: Not that much shorter.
  %           error(irf.assert.ASSERTION_EMID, msg)
  %       ==> irf.assert.error(msg)
  %           assert(... , irf.assert.ASSERTION_EMID, msg)
  %       ==> irf.assert.assert(... , msg)
  %   CON: Gets method ~irf.assert.assert() is a bad name.
  %       CON: Private function so matters less.
  %   CON: Potentially slower?
  %       TODO-NI
  %       CON: error only applies to an internal assert function. Internal "error" function only triggered when assertion
  %            already has failed.
  %   PRO: Easier to be consistent.
  %   PRO: Clearer.
  %       PRO: Fewer arguments.
  %
  %
  %
  % PROPOSAL: Add argument for name of variable (used as argument) so that can automatically print better error messages.
  % PROPOSAL: Optional last arguments for (1) error message, or (2) error msg ID + error message.
  %   CON: Can conflict with other string arguments.
  %       CON-PROPOSAL: For assertions with arbitrary number of arguments
  %       (varargin), use cell array instead.
  %           Ex: sizes({xSc, [-1, 1], ySc, [-1, 1]}, 'xpilot.BadSpacecraftShape', 'Illegal spacecraft shape specification.')
  %           PRO: Only sizes() uses varargin 2020-10-07.
  % PROPOSAL: Same syntax as for last arguments of MATLAB function assert(), i.e.
  %   (1) error message
  %   (2) error message ID + error message (sprintf pattern) + (optional) variable values
  %       for error message (sprintf pattern)
  %   PROPOSAL: Internal assertion error function can process varargin for every
  %             assertion function.
  %
  %
  % PROPOSAL: Have way of obtaining whether assertion is satisfied, without throwing error.
  %   PRO: Can use assertion methods for raising customized exceptions (not just assertion exceptions), e.g. for UI
  %        errors.
  %   PRO: Useful for creating more compact custom-made assertions.
  %       Ex: assert(isscalar(padValue) || is_castring(padValue))
  %       Ex: assert(<castring> || <struct>)
  %       Ex: Checking settings values.
  %           Ex: assert(ischar(defaultValue) || isnumeric(defaultValue) || is_castring_set(defaultValue))
  %   PRO: Can construct custom-made assertion that can not be created through combinations of pre-existing assertion
  %       functions.
  %       Ex: assert(isscalar(padValue) || is_castring(padValue))
  %   PRO: Can use assertion methods for checking state/conditions (to execute code; not raise errors).
  %       Ex: if <castring> elseif <castring_set> else ... end
  %
  %   CON: The dividing line between (a) assertion (which raise error) on the inputs of the condition function, and (b)
  %        the result of the function itself becomes blurry.
  %        IMPORTANT: Assertion functions should have no effect (barring bugs in the outside code). They can assert more
  %        and less as long as they assert a subset of what can be asserted. If an assertion is tightened or loosened that
  %        way, then that (ideally) has no effect. For a condition function, tightening/loosening the definition of a
  %        condition may have consequences.
  %       Ex: Confusing arguments for value to be tested, and parameters to the test/assertion.
  %       Ex: ~is_castring_set: Should the function return false, or raise assertion error on input if value is
  %           (1) not a cell array, or
  %           (2) cell array of non-strings, or
  %           (3) cell array of non-unique strings ?
  %       Ex: is_castring_in_set(s, strSet): Should the function return false when
  %           (1) s is string in strSet, when strSet is a cell array of non-unique strings, or
  %           (2) s is not a string ?
  %       Ex: ~have_sizes(): sizeConstraints is not a numeric 1D array.
  %       Ex: ~castring_sets_equal()
  %       Ex: ~file_exists(path)
  %           False if path is not a string?
  %       PROPOSAL: Other naming convention that makes it clear what is assumed to be true, and what is properly asserted
  %           Ex: castring_set --> ca_of_strings_is_set
  %
  %   PROPOSAL: Condition results should correspond to the exact condition implied in the function name, and
  %             OPTIONALLY use internal assertions for other conditions implied in the name.
  %       Ex: ~is_castring_in_set(s, strSet)
  %           Optionally assert that s is castring, and strSet is set. Return result, assuming those are true.
  %           CON: Sounds wrong.
  %       Ex: castring_sets_equal(set1, set2)
  %           Optionally assert that set1 and set2 are sets.
  %       PROPOSAL: Interpret condition from the adjective or verb.
  %       PROPOSAL: Assume that arguments are what they "should be", in isolation.
  %
  %   PROBLEM: How handle assertions functions that also return values.
  %       Ex: sizes()
  %           PROPOSAL: Condition function identical. Return boolean ahead of all other return values.
  %               Ex: [satisfied, nRecords] = ~have_sizes(zVar1, [-1], zVar2, [-1, 2048])
  %
  %   PROPOSAL: Have assertion methods return a true/false value for assertion result.
  %             If a value is returned, then raise no error.
  %       CON: sizes() already has return values (and aritrarily many).
  %       PRO: Can be used with regular assert statements.
  %           PRO: MATLAB's assert can be used for raising exception with customized error message.
  %
  %   PROPOSAL: Have analogous but SEPARATE "CONDITION FUNCTIONS".
  %       PROPOSAL: Naming scheme  is_*, has_*, have_*
  %           CON: Not perfect since also sounds like assertions.
  %               Ex: assert.have_sizes, assert.is_castring
  %           CON: Does not work for all assertion functions
  %               Ex: castring_sets_equal --> are_castring_sets_equal? castring_sets_are_equal?
  %       CON: "assert" is a bad name for such a class.
  %       PROPOSAL: Create class with collection of standardized non-trivial
  %                 "condition functions", used by this "assert" class.
  %                 Use an analogous naming scheme.
  %           CON: From the point of view of other code, there is nothing special about the "conditions functions".
  %           	 They are just regular functions. Should not be specially grouped.
  %           PROPOSAL: Class name "cond".
  %       PROPOSAL: Only do so for well-defined conditions.
  %           Ex: Not castring.
  %
  %
  %
  % PROPOSAL: Functions for asserting line breaks.
  %   TODO-DEC: Which set of functions.
  %       PROPOSAL: Assert ending LF (assert not ending CR+LF).
  %       PROPOSAL: Assert all linebreaks are LF (no CR+LF).
  %       PROPOSAL: Assert all linebreaks are LF (no CR+LF). Require ending linebreak.
  %
  % PROPOSAL: Assertion functions for MATLAB's date vectors.
  %   NOTE: Variants with 3 and 6 components
  %       PROPOSAL: datevec3, datevec6
  %       PROPOSAL: datevec(dv, nComp)    % nComp = 3,6
  %   NOTE: Variants with 1 or many rows.
  %       PROPOSAL: datevec(dv, nComp, oneManyRows)
  %
  % PROPOSAL: Assertion like path_is_available() but for string pattern.
  %   Ex: Assertion that file (dataset? document?) of any version does not exist
  %       at location.
  %       Ex: Dataset does not exist in /data/solo/soar/.
  %


  properties(Constant, Access=private)
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
        error(irf.assert.ASSERTION_EMID, ...
          'Expected castring (0x0, 1xN char array) is not char.')

      elseif ~(isempty(s) || size(s, 1) == 1)
        error(irf.assert.ASSERTION_EMID, ...
          'Expected castring (0x0, 1xN char array) has illegal dimensions.')
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
    %
    function castring_regexp(s, regexp)
      assert(ischar(s))
      if ~any(irf.str.regexpf(s, regexp))
        error(irf.assert.ASSERTION_EMID, ...
          ['String "%s" (in its entirety) does not match any of', ...
          ' the specified regular expressions.'], s)
      end
    end



    % Argument is a cell matrix of UNIQUE strings.
    function castring_set(s)
      % NOTE: Misleading name, since does not check for strings.

      assert(iscell(s), ...
        irf.assert.ASSERTION_EMID, ...
        'Expected cell array of unique strings, but is not cell array.')

      % IMPLEMENTATION NOTE: For cell arrays, "unique" requires the
      % components to be strings. Therefore does not check (again), since
      % it is probably slow.
      assert(numel(unique(s)) == numel(s), ...
        irf.assert.ASSERTION_EMID, ...
        'Expected cell array of unique strings, but not all strings are unique.')
    end



    % NOTE: Does not care about dimensions.
    function castring_sets_equal(set1, set2)
      % NOTE/BUG: Does not require sets to have internally unique strings.

      % NOTE: setxor() requires cell arrays of strings(?).
      if ~isempty(setxor(set1, set2))
        error(irf.assert.ASSERTION_EMID, ...
          'The two string sets are not equivalent.')
      end
    end



    function castring_in_set(s, strSet)
      % PROPOSAL: Abolish
      %   PRO: Unnecessary since can use assert(ismember(s, strSet)).
      %       CON: This gives better error messages for string not being string,
      %            for string set not being string set.

      irf.assert.castring_set(strSet)
      irf.assert.castring(s)

      if ~ismember(s, strSet)
        error(irf.assert.ASSERTION_EMID, ...
          'Expected string in string set is not in set.')
      end
    end



    % NOTE: Can also be used for checking supersets.
    % NOTE: Both string sets and numeric sets
    function subset(strSubset, strSuperset)

      if ~irf.utils.subset(strSubset, strSuperset)
        error(irf.assert.ASSERTION_EMID, ...
          'Expected subset/superset is not a subset/superset.')
      end
    end



    % Check that numeric/logical array contains only unique values.
    % NaN, Inf, -Inf count as a unique values.
    function number_set(v)
      % NOTE: number_set analogous to castring_set.
      % PROPOSAL: Better name considering the accepted MATLAB classes.

      % IMPLEMENTATION NOTE: Special cases for "unique".
      %   NaN:      Every NaN counts as unique (i.e. different from other NaN).
      %   In, -Inf: Inf counts as equal to itself, -Inf counts as equal to
      %              itself.
      % Must therefore count the number of NaNs.
      % NOTE: Works for logical, (individual) character arrays, but not for cell arrays of strings.
      assert(sum(isnan(v)) <= 1, ...
        irf.assert.ASSERTION_EMID, ...
        'Array does not contain only unique numbers. It contains multiple NaN.')

      % IMPLEMENTATION NOTE: Also works for strings (but the NaN check above does not).
      assert(numel(unique(v)) == numel(v), ...
        irf.assert.ASSERTION_EMID, ...
        'Array does not contain only unique numbers.')
    end



    function scalar(x)
      if ~isscalar(x)
        error(irf.assert.ASSERTION_EMID, ...
          'Variable is not scalar as expected.')
      end
    end



    % Either regular file or symlink to regular file (i.e. not directory or
    % symlink to directory).
    % Cf. path_is_available().
    function file_exists(filePath)
      if ~(exist(filePath, 'file') == 2)
        error(irf.assert.ASSERTION_EMID, ...
          ['Expected existing regular file (or symlink to regular', ...
          ' file) "%s" can not be found.'], filePath)
      end
    end



    function dir_exists(dirPath)
      if ~exist(dirPath, 'dir')
        error(irf.assert.ASSERTION_EMID, ...
          'Expected existing directory "%s" can not be found.', ...
          dirPath)
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
        error(irf.assert.ASSERTION_EMID, ...
          ['Path "%s" which was expected to point to nothing,', ...
          ' actually points to a file or directory.'], path)
      end
    end



    % ARGUMENTS
    % ==========
    % requiredFnSet
    %       Cell array of required field names.
    % optionalFnSet
    %       One of two alternatives.
    %       (a) Cell array of optional field names (i.e. allowed, but not
    %           required)
    %       (b) String constant 'all'
    %           All fieldnames are allowed but not required. This is likely
    %           only meaningful when requiredFnSet is non-empty (not a
    %           requirement).
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
      %   TODO-DEC: How specify fieldnames? Can not use cell arrays recursively.
      %   PROPOSAL: Define other, separate assertion method.
      %   PROPOSAL: Tolerate/ignore that structs are array structs.
      %   PROPOSAL: struct(S, {'PointA.x', 'PointA.y'}, {'PointA.z'})
      %   PROPOSAL: struct(S, {'PointA', {'x', 'y'}}, {'PointA', {'z'}})   % Recursively
      %       Cell array means previous argument was the parent struct.
      %   PROPOSAL: struct(S, {'name', 'ReqPointA', {{'reqX', 'reqY'}, {'optZ'}}},
      %           {'OptPointB', {{'reqX', 'reqY'}, {'optZ'}}})
      %       Cell array means previous argument was the parent struct.
      %       Groups together required and optional with every parent struct.
      %       PRO: Optional fields can be structs with both required and optional fields, recursively.
      %       PRO: Can be implemented recursively(?).
      %   TODO-NI: Required & optional is well-defined?
      %   CON: Rarely needed.
      %       CON: Maybe not
      %           Ex: Settings structs
      %   CON-PROPOSAL: Can manually call irf.assert.struct multiple times, once for each substruct,
      %                 instead (if only required field names).
      %
      % PROPOSAL: Assertion: Intersection requiredFnSet-optionalFnSet is empty.
      %
      % PROPOSAL: Separate out condition function.
      %   PRO: Useful for distinguishing between multiple "struct formats".

      assert(isstruct(S))
      structFnSet          = fieldnames(S);

      missingRequiredFnSet = setdiff(requiredFnSet, structFnSet);

      % disallowedFnSet = ...
      if iscell(optionalFnSet)
        disallowedFnSet = setdiff(...
          setdiff(structFnSet, requiredFnSet), ...
          optionalFnSet);
      elseif isequal(optionalFnSet, 'all')
        disallowedFnSet = {};
      else
        error(irf.assert.ASSERTION_EMID, ...
          ['Illegal optionalFnSet argument.', ...
          ' Is neither cell array or string constant "all".'])
      end

      % Give error, with an actually useful error message.
      if ~isempty(missingRequiredFnSet) || ~isempty(disallowedFnSet)
        missingRequiredFnListStr = strjoin(missingRequiredFnSet, ', ');
        disallowedFnListStr      = strjoin(disallowedFnSet,      ', ');

        error(irf.assert.ASSERTION_EMID, ...
          ['Expected struct has the wrong set of fields.', ...
          '\n    Missing fields:           %s', ...
          '\n    Extra (forbidden) fields: %s'], ...
          missingRequiredFnListStr, disallowedFnListStr)
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
    % NOTE: Can not handle function handle that does not point to existing
    % function(!)
    %
    % NOTE: Method's USEFULNESS IS LIMITED due to below reasons:
    %   NOTE: Can not distinguish between functions that
    %       (1) only accept exactly N arguments, and
    %       (2) accept a variable-number of arguments, including N arguments
    %           (using varargin).
    %       Analogous problem with nargout.
    %   NOTE: Empirically, nargout(anonymousFunction) == -1, always.
    %
    function func(funcHandle, nArgin, nArgout)
      % PROPOSAL: Make nArgin/nArgout able to simultaneously accept multiple nargin/nargout values.
      %   Ex: Accept nargin  = 3, -1,-2,-3 (accept all functions that seemingly can accept three arguments)
      %   Ex: Accept nargout = 3, -1,-2,-3 (accept all functions that seemingly can return three values)
      %   PROPOSAL: Submit exact set (numeric array) of accepted values.
      % PROPOSAL: Do not submit/assert any specific number of arguments or
      %           return values.
      %   PROPOSAL: Separate assertion function.

      if ~isa(funcHandle, 'function_handle')
        error(irf.assert.ASSERTION_EMID, ...
          'Expected function handle is not a function handle.')
      end
      if nargin(funcHandle) ~= nArgin
        error(irf.assert.ASSERTION_EMID, ...
          ['Expected function handle ("%s") has the wrong number', ...
          ' of input arguments. nargin()=%i, nArgin=%i'], ...
          func2str(funcHandle), nargin(funcHandle), nArgin)
      elseif nargout(funcHandle) ~= nArgout
        % NOTE: MATLAB actually uses term "output arguments".
        error(irf.assert.ASSERTION_EMID, ...
          ['Expected function handle ("%s") has the wrong number of output', ...
          ' arguments (return values). nargout()=%i, nArgout=%i'], ...
          func2str(funcHandle), nargout(funcHandle), nArgout)
      end
    end



    % Assert v has a non-one size in at most one dimension.
    % NOTE: Excludes 0x0, e.g. [], {}, and ''.
    %
    % NOTE: MATLAB's "isvector" function uses different criterion which
    % excludes length in third or higher dimension.
    %   isvector(ones(0,1)) == 1
    %   isvector(ones(1,0)) == 1
    %   isvector(ones(1,1,3)) == 0
    %   isvector(ones(1,1,0)) == 0
    %
    function vector(v)
      % PROPOSAL: Optional extra argument that specifies the length.
      %   CON: Functionality should be incorporated into irf.assert.sizes() (?).
      % PROPOSAL: Better name.
      %   PRO: "vector" is ambiguous w.r.t. number of dimensions.
      %   PROPOSAL: Name that implies 1D vector: vector_1D, vec_1D
      %   PROPOSAL: Name that implies exactly one non-zero size dimension.
      %       PROPOSAL: true_vector, strict_vector, strict_vec_1D, vec_strict_1D
      % PROPOSAL: Configurable to permit []={}=0x0 specifically
      %   PROPOSAL: Separate policy constant to permit. '0x0', 'permit 0x0'.
      %   PROPOSAL: Separate method to permit.
      %       PROPOSAL: common_vec_1D, nonstrict_vec_1D.

      if sum(size(v) ~= 1) > 1
        % CASE: More than one dimension has size > 1.
        sizeStr = sprintf('%ix', size(v));
        sizeStr = sizeStr(1:end-1);   % ~Hack

        error(irf.assert.ASSERTION_EMID, ...
          'Expected vector, but found variable of size %s.', sizeStr)
      end
    end



    % Assert the sizes of one or multiple variables.
    % See irf.ds.sizes for arguments and return values.
    %
    % NOTE: One (sometimes) needs to add semicolon to end of row, since has
    %       return values.
    % NOTE: Can use the return values to assert conditions like n>10 (after
    %       call to this function), which this function itself can not.
    % NOTE: Returning useful values in principle make this function have
    %       both assertion and non-assertion functionality.
    %   Ex: Returning number of CDF records while simultaneously asserting
    %       consistent sizes of multiple (zv/MATLAB) variables.
    % NOTE: Can be used to verify that sizes of variables are equal, without
    %       constraints on numebr of dimensions.
    %           Ex: irf.assert.sizes(A, size(B))
    %       NOTE: This takes care of problem of trailing ones (must be
    %             normalized to be equal for equal sizes, e.g. removed) when
    %             comparing size(A) and size(B).
    %
    function [varargout] = sizes(varargin)
      % See irf.ds.sizes() for BOGIQ.

      varargout = cell(1, nargout);
      [condSatisfied, varargout{:}] = irf.ds.sizes(varargin{:});

      assert(condSatisfied, irf.assert.ASSERTION_EMID, ...
        'Variable sizes do not satisfy specified constraints.')
    end



    % Assert that all values in a matrix are identical.
    %
    % NOTE: Cf isequal() with n>=2 arguments. Does not generalize to n=1,
    % n=0, or higher-dimensional matrices.
    % argument.
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
    %     NOTE: Empty or scalar matrices are accepted.
    %
    function all_equal(v)
      % TODO-DEC: How handle NaN?
      %   PROPOSAL: Count NaN as equal to itself.

      nUniques = numel(unique(v(:)));    % NOTE: Make 1D vector.
      nTotal   = numel(v);

      if (nUniques ~= 1) && (nTotal >= 1)
        error(irf.assert.ASSERTION_EMID, ...
          ['Expected vector of identical values, but found %i', ...
          ' unique values out of a total of %i values.'], ...
          nUniques, nTotal)
      end
    end



  end    % methods(Static)



end    % classdef
