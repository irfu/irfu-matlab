%
% Class for static methods for creating assertions.
%
% NOTE: MATLAB already (at least MATLAB R2009a) has a function "assert" which is useful for simpler cases.
%
%
% POLICY
% ======
% Functions should be named as propositions (when including the class name "assert") which are true if the assertion function does not yield error.
% "castring" refers to arrays of char, not the concept of "strings" which begins with MATLAB 2017a and later.
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
%     vs EJ_library.utils.assertions.castring_in_set(datasetType, {'DERIV1'})
% Ex: assert(any(strcmp(s, {'DERIV1', 'EDDER'})))
%     vs EJ_library.utils.assertions.castring_in_set(datasetType, {'EDDER', 'DERIV1'})
% Ex: assert(isstruct(s) && isempty(setxor(fieldnames(s), {'a', 'b', 'c'})))
%     vs EJ_library.utils.assertions.is_struct_w_fields(s, {'a', 'b', 'c'})
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
% TODO-DECISION: Use assertions on (assertion function) arguments internally?
% PROPOSAL: Add argument for name of argument so that can print better error messages.
% PROPOSAL: Optional error message (string) as last argument to ~every method.
%   CON: Can conflict with other string arguments.
%       Ex: Method "struct".
% PROPOSAL: Optional error message identifier as second-last argument to ~every method (see "error").
%   CON: Can conflict with other string arguments.
%
% PROPOSAL: Assertion for checking that multiple variables have same size in specific dimensions/indices.
%   PROPOSAL: Same dimensions in all dimensions except those specified.
%       PRO: Better handling of "high dimensions to infinity".
%   PROPOSAL: Check on all fields in struct.
%       Ex: SSL (+KVPL?)
%
% PROPOSAL: Redefine class as collection of standardized non-trivial "condition functions", used by an "assert" class.
% PROPOSAL: Have methods return a true/false value for assertion result. If a value is returned, then raise no assertion error.
%   PRO: Can be used with regular assert statements. (Not entirely sure why this would be useful.)
%       PRO: MATLAB's assert can be used for raising exception with customized error message.
%   PRO: Can use assertion methods for checking state/conditions (without raising errors).
%       Ex: if <castring> elseif <castring_set> else ... end
%   PRO: Can use assertion methods for raising customized exceptions (not just assertion exceptions), e.g. for UI
%        errors.
%   PRO: Can combine multiple assertion conditions into one assertion.
%       Ex: assert(<castring> || <struct>)
%
% PROPOSAL: Static variable for error message identifier.
%   PRO: Can set differently in BICAS.
%
% PROPOSAL: Move out of "utils" package.
%   PRO: Is used very, very often.



    properties(Constant)
        ERROR_MSG_ID = 'assert:Assertion'
    end

    
    
    methods(Static)
        
        % NOTE: Empty string literal '' is 0x0.
        % NOTE: Accepts all empty char arrays, e.g. size 2x0 (2 rows).
        function castring(s)
            % PROPOSAL: Only accept empty char arrays of size 1x0 or 0x0.
            if ~ischar(s)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected castring (0x0, 1xN char array) is not char.')
            elseif ~(isempty(s) || size(s, 1) == 1)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected castring (0x0, 1xN char array) has wrong dimensions.')
            end
        end
        
        
        
        % Assert that ENTIRE string matches a regexp or any in a cell array of regular expressions.
        % NOTE: If regex is an empty cell array, then assertion fails.
        % NOTE: Will permit empty strings to match a regular expression.
        function castring_regexp(s, regexp)
            if ~any(EJ_library.utils.regexpf(s, regexp))
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'String "%s" (in its entirety) does not match any of the specified regular expressions.', s)
            end
        end
        
        
        
        % Cell matrix of UNIQUE strings.
        function castring_set(s)
            % NOTE: Misleading name, since does not check for strings.
            
            if ~iscell(s)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected cell array of unique strings, but is not cell array.')
                
            % IMPLEMENTATION NOTE: For cell arrays, "unique" requires the components to be strings. Therefor does not
            % check (again), since probably slow.
            elseif numel(unique(s)) ~= numel(s)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected cell array of unique strings, but not all strings are unique.')
            end
        end

        
        
        function castring_in_set(s, strSet)
        % PROPOSAL: Abolish
        %   PRO: Unnecessary since can use assert(ismember(s, strSet)).
        %       CON: This gives better error messages for string not being string, for string set not being string set.
            import EJ_library.*
        
            utils.assert.castring_set(strSet)
            utils.assert.castring(s)
            
            if ~ismember(s, strSet)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected string in string set is not in set.')
            end
        end
        
        
        
        % NOTE: Can also be used for checking supersets.
        function subset(strSubset, strSet)
            % PROPOSAL: Name without "string", since does not check for strings.
            %   PROPOSAL: Change name?
            
            % NOTE: all({}) == all([]) == true
            if ~all(ismember(strSubset, strSet))
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected subset is not a subset.')
            end
        end
        
        
        
        function scalar(x)
            if ~isscalar(x)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Variable is not scalar as expected.')
            end
        end
        
        
        
        % Either regular file or symlink to regular file (i.e. not directory or symlink to directory).
        function file_exists(filePath)
            if ~(exist(filePath, 'file') == 2)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected existing regular file (or symlink to regular file) "%s" can not be found.', filePath)
            end
        end
        
        
        
        function dir_exists(dirPath)
            if ~exist(dirPath, 'dir')
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected existing directory "%s" can not be found.', dirPath)
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
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Path "%s" which was expected to point to nothing, actually points to a file/directory.', path)
            end
        end
        
        
        
        % Struct with certain set of fields.
        % 
        % ARGUMENTS
        % =========
        % varargin :
        %   <Empty>    : Require exactly                   the specified set of fields.
        %   'subset'   : Require subset   of (or equal to) the specified set of fields.
        %   'superset' : Require superset of (or equal to) the specified set of fields.
        %
        % NOTE: Does NOT assume 1x1 struct. Can be matrix.
%         function struct(S, fieldNamesSet, varargin)
%             % NEED: Be able to specify
%             %   struct with exact set of fieldnames
%             %   struct with some required fieldnames, and some optional fieldnames
%             %   struct with only optional fieldnames                                               <=> fieldnames are a SUBSET of specified set.
%             %   struct with with some required fieldnames, and all other fieldnames being optional <=> fieldnames are a SUPERSET of specified fieldnames
%             % 
%             % PROPOSAL: Replace superset, subset with clearer keywords.
%             %   PROPOSAL: require, permit
%             % PROPOSAL: Change arguments for fieldnames. Always two.
%             %   PRO: More robuts w.r.t. backward-compatibility in the case of future changes(?)
%             %   PROPOSAL: requiredFnSet, optionalFnSet
%             %       PROPOSAL: optionalFnSet can have special value ('all') to mean that all fieldnames are optional.
%             %
%             % PROPOSAL: Delete.
%             
%             warning('Using deprecated function EJ_library.utils.assert.struct.')
%             
%             import EJ_library.*
%             
%             if isempty(varargin)   %numel(varargin) == 1 && isempty(varargin{1})
%                 checkType = 'exact';
%             elseif numel(varargin) == 1 && strcmp(varargin{1}, 'subset')
%                 checkType = 'subset';
%             elseif numel(varargin) == 1 && strcmp(varargin{1}, 'superset')
%                 checkType = 'superset';
%             else
%                 error(EJ_library.utils.assert.ERROR_MSG_ID, 'Illegal argument')
%             end
%             
%             if ~isstruct(S)
%                 error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected struct is not struct.')
%             end
%             utils.assert.castring_set(fieldNamesSet)    % Abolish?
%             
%             missingFnList = setdiff(fieldNamesSet, fieldnames(S));
%             extraFnList   = setdiff(fieldnames(S), fieldNamesSet);
%             
%             switch(checkType)
%                 case 'exact'
%                     if (~isempty(missingFnList) || ~isempty(extraFnList))
%                         
%                         missingFnListStr = strjoin(missingFnList, ', ');
%                         extraFnListStr   = strjoin(extraFnList,   ', ');
%                         
%                         error(EJ_library.utils.assert.ERROR_MSG_ID, ['Expected struct has the wrong set of fields.', ...
%                             '\n    Missing fields:           %s', ...
%                             '\n    Extra (forbidden) fields: %s'], missingFnListStr, extraFnListStr)
%                     end
%                 case 'subset'
%                     if ~isempty(extraFnList)
%                         
%                         extraFnListStr   = strjoin(extraFnList,   ', ');
%                         error(EJ_library.utils.assert.ERROR_MSG_ID, ['Expected struct has the wrong set of fields.', ...
%                             '\n    Extra (forbidden) fields: %s'], extraFnListStr)
%                     end
%                 case 'superset'
%                     if ~isempty(missingFnList)                        
%                         missingFnListStr = strjoin(missingFnList, ', ');
%                         error(EJ_library.utils.assert.ERROR_MSG_ID, ['Expected struct has the wrong set of fields.', ...
%                             '\n    Missing fields:           %s'], missingFnListStr)
%                     end
%             end
%         end



        % Replacement for "EJ_library.utils.assert.struct".
        %
        % ARGUMENTS
        % ==========
        % requiredFnSet : Cell array of required field names.
        % optionalFnSet : (a) Cell array of optional field names (i.e. allowed, but not required)
        %                 (b) String constant 'all' : All fieldnames are allowed but not required. This is likely only
        %                     meaningful when requiredFnSet is non-empty (not a requirement).
        %
        function struct2(S, requiredFnSet, optionalFnSet)
            % PROPOSAL: Have it apply to a set of strings (e.g. fieldnames), not a struct as such.
            % PROPOSAL: Let optionalFnSet be optional (empty by default).
            %   PRO: Shorter for the most common case.
            %   PRO: Backward-compatibility with some of the syntax for struct(predecessor assertion function).
            %   CON: Bad for future extensions of function.
            %
            % PROPOSAL: Recursive structs field names.
            %   TODO-DECISION: How specify fieldnames? Can not use cell arrays recursively.
            %   PROPOSAL: Define other, separate assertion method.
            %   CON: Rarely needed.
            %   CON-PROPOSAL: Can manually call EJ_library.utils.assert.struct2 multiple times, once for each substruct,
            %                 instead (if only required field names).
            % PROPOSAL: Rename struct2-->struct.
            % PROPOSAL: Assertion: Intersection requiredFnSet-optionalFnSet is empty.
            
            structFnSet          = fieldnames(S);
            
            missingRequiredFnSet = setdiff(requiredFnSet, structFnSet);
            
            % disallowedFnSet = ...
            if iscell(optionalFnSet)
                disallowedFnSet = setdiff(setdiff(structFnSet, requiredFnSet), optionalFnSet);
            elseif isequal(optionalFnSet, 'all')
                disallowedFnSet = {};
            else
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Illegal optionalFnSet argument. Is neither cell array or string constant "all".')
            end
            
            % Give error, with an actually useful error message.
            if ~isempty(missingRequiredFnSet) || ~isempty(disallowedFnSet)
                missingRequiredFnListStr = strjoin(missingRequiredFnSet, ', ');
                disallowedFnListStr      = strjoin(disallowedFnSet,      ', ');

                error(EJ_library.utils.assert.ERROR_MSG_ID, ['Expected struct has the wrong set of fields.', ...
                    '\n    Missing fields:           %s', ...
                    '\n    Extra (forbidden) fields: %s'], missingRequiredFnListStr, disallowedFnListStr)
            end
        end



        % NOTE: Can not be used for an assertion that treats functions with/without varargin/varargout.
        %   Ex: Assertion for functions which can ACCEPT (not require exactly) 5 arguments, i.e. incl. functions which
        %       take >5 arguments.
        % NOTE: Not sure how nargin/nargout work for anonymous functions. Always -1?
        % NOTE: Can not handle: is function handle, but does not point to existing function(!)
        function func(funcHandle, nArgin, nArgout)
            if ~isa(funcHandle, 'function_handle')
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected function handle is not a function handle.')
            end
            if nargin(funcHandle) ~= nArgin
                error(EJ_library.utils.assert.ERROR_MSG_ID, ...
                    'Expected function handle ("%s") has the wrong number of input arguments. nargin()=%i, nArgin=%i', ...
                    func2str(funcHandle), nargin(funcHandle), nArgin)
            elseif nargout(funcHandle) ~= nArgout
                % NOTE: MATLAB actually uses term "output arguments".
                error(EJ_library.utils.assert.ERROR_MSG_ID, ...
                    'Expected function handle ("%s") has the wrong number of output arguments (return values). nargout()=%i, nArgout=%i', ...
                    func2str(funcHandle), nargout(funcHandle), nArgout)
            end
        end
        
        
        
        function isa(v, className)
            if ~isa(v, className)
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected class=%s but found class=%s.', className, class(v))
            end
        end
        
        

        % Assert v has a non-one size in at most one dimension.
        %
        % NOTE: MATLAB's "isvector" function uses different criterion which excludes numel(v) == 0, and length in third
        % or higher dimension.
        function vector(v)
            dims = unique(size(v));
            dims(dims==1) = [];
            if numel(dims) > 1
                sizeStr = sprintf('%ix', size(v));
                error(EJ_library.utils.assert.ERROR_MSG_ID, 'Expected vector, but found variable of size %s.', sizeStr(1:end-1))
            end
        end
        
        
        
        % Assert that all values in a matrix are identical. Useful for e.g. checking that sizes of vectors are
        % identical.
        %
        % NOTE: Empty matrices are accepted.
        %
        % Works on:
        % - matrix of numbers
        % - matrix of characters
        % - cell array of strings
        %
        % Does not work on:
        % - cell array of numbers
        % 
        function all_equal(v)
           nUniques = numel(unique(v(:)));    % NOTE: Make 1D vector.
           nTotal   = numel(v);
           if (nUniques ~= 1) && (nTotal >= 1)
               error(EJ_library.utils.assert.ERROR_MSG_ID, ...
                   'Expected vector of identical values, but found %i unique values out of a total of %i values.', ...
                   nUniques, nTotal)
           end
        end
        
    end    % methods
end    % classdef
