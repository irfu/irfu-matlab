%
% Class which defines a type of automatic test: Check whether an arbitrary
% function (to be tested) produces the specified expected return result
% (deterministic).
%
% NOTE: The exact function comparing expected with actual return result can be
% specified, and thus the comparison can be relaxed, e.g. to
%   * allow for numeric imprecision
%   * only check the results of some return values
%
% NOTE: The default comparison function is "isequaln". This is not ideal since
% it does not care about matlab class, e.g.
% isequal(0, false) == isequal('A', 65) == true.
%
%
% TERMINOLOGY, ABBREVIATIONS
% ==========================
% expected = exp : 
% actual   = act : 
%       Something that the function was EXPECTED TO/ACTUALLY do/did or returned
%       when called.
% input          : Arguments to function (all combined).
% output         : Function return values (all combined).
%
%
% NOTE: HOW TO TEST METHODS
% =========================
% Use function handles to the methods as below:
%   Class instance (not static?) methods (not constructors):
%       @methodNameExcludingClassOrPackages, where the first input argument is
%       the object/class instance.
%   Class instance constructor
%       @classNameIncludingPackages
%
%
% NOTE: For the case of expected exception, the test is not able to make the
% call with a specified number of return values (output).
%
% BUG: Expecting too few return values (having too many actual return values) is
% mistakenly accepted.
%
%
% Initially created 2019-01-21 by Erik P G Johansson.
%
classdef CompareFuncResult
    %
    % PROPOSAL: Always specify number of return values, also for expected exception.
    %   PRO: nargout could influence execution.
    %   CON: Would require one more argument (constructor) for expected exception.
    %   PROPOSAL: Expected exception (only) constructor: func, input, expOutputOrNOutput, exceptionName
    %   PROPOSAL: Expected exception (only) constructor: func, input, expOutputOrExcName, nOutput
    %   PROPOSAL: Expected exception        constructor: func. input, nOutput, expOutput
    %             Expected non-exception    constructor: func. input, nOutput, expOutputOrExcName   % Assertion: nOutput=...
    % PROPOSAL: First constructor argument = string ==> Use as name for test.
    %   PRO: Useful for identifying test which fails, and selecting that only.
    %
    % PROPOSAL: Specify which return value(s) differ.
    %
    % PROPOSAL Reorganize as testing a function call.
    %   Does arbitrary consistency tests on the combination of input, actual output,
    %   and test-supplied data. Actual output equalling test-supplied data
    %   (=expected output) is only one special case.
    %   NOTE: Cf python counterpart.
    %   PRO: All the handling of exceptions (expeected vs actual, exception class) is shared between all tests of
    %        functions.
    %   CON: The format of the log messages is dependent on the type of test.
    %       NOTE: Only the automatic search for difference between
    %             actual output and user-supplied data/expected output.
    %   PROPOSAL: Implement this class as subclass to superclass which is
    %             generic according to proposal. This class would only implement
    %             comparisons between expected & actual output.
    %
    properties(SetAccess=private,GetAccess=public)
        func
        input
        expOutput
        expExceptionName
        Settings
    end



    methods(Access=public)
        
        % Constructor
        %
        %
        % ARGUMENTS
        % =========
        % func               
        %       Function handle
        % input  
        %       Cell array of function arguments
        % expOutputOrExcName
        %       Either
        %           (1) cell array of expected return values, or
        %           (2) character string naming the expected Exception class.
        %       NOTE: The default exception class is "MException".
        % varargin
        %       Settings as interpreted by
        %       EJ_library.utils.interpret_settings_args.
        %       'equalsFunc' : Function handle
        %           areEqual = equalsFunc(expOutputCa, actOutputCa)
        %       'keyboardOnUnequal' : 0/1, boolean.
        %           Whether to stop execution with "keyboard" so the user can
        %           manually compare and inspect the expected and actual output.
        %
        function obj = CompareFuncResult(func, input, expOutputOrExcName, varargin)
            % ASSERTION
            assert(iscell(input), 'Argument "input" is not a cell array')
            
            obj.func  = func;
            obj.input = input;
            if iscell(expOutputOrExcName)
                obj.expOutput        = expOutputOrExcName;
                obj.expExceptionName = [];
            elseif ischar(expOutputOrExcName)
                obj.expOutput        = [];
                obj.expExceptionName = expOutputOrExcName;
            else
                % ASSERTION
                error(...
                    ['expOutputOrExcName is neither', ...
                    ' (1) cell array, nor (2) string.'])
            end
            
            % MATLAB R2016a: "isequalwithequalnans is not recommended. Use
            % ISEQUALN instead."
            DEFAULT_SETTINGS.equalsFunc        = @isequaln;
            DEFAULT_SETTINGS.keyboardOnUnequal = 0;
            
            obj.Settings   = EJ_library.utils.interpret_settings_args(...
                DEFAULT_SETTINGS, varargin);
            EJ_library.assert.struct(...
                obj.Settings, fieldnames(DEFAULT_SETTINGS), {})
        end



        function TestData = run(obj)
            try
                %---------------------------------------------------------------
                % IMPLEMENTATION NOTE: It is non-trivial to make this code work
                % for all combinations of (1) return values (none, not none,
                % variable number), and (2) exception, non-exception.
                %---------------------------------------------------------------
                if isempty(obj.expOutput) || ~isempty(obj.expExceptionName)
                    % CASE: Expecting (1) no output, or (2) exception.
                    
                    %===========
                    % Call test
                    %===========
                    obj.func(obj.input{:});
                    % Could be incorrect, but only for the case of expected
                    % exception, and actual non-exception
                    actOutput = {};
                else
                    actOutput = cell(size(obj.expOutput));
                    %==========================================================
                    % Call test
                    % ---------
                    % IMPLEMENTATION NOTE: Important with
                    % (1) square brackets around return result (error for zero
                    %     return values)
                    % (2) variable which captures return values is pre-defined
                    %     with correct size.
                    %==========================================================
                    [actOutput{:}] = obj.func(obj.input{:});
                end
                actException = [];
            catch Exc
                actOutput    = [];
                actException = Exc;
            end

            if ~isempty(actException)
                %========================
                % CASE: Actual exception
                %========================
                
                if ~isempty(obj.expExceptionName)
                    %==========================
                    % CASE: Expected exception
                    %==========================
                    if strcmp(class(actException), obj.expExceptionName)
                        %=================================================
                        % CASE: Actual and expected exception classes are
                        %       identical.
                        %=================================================
                        TestData.resultDescrText = 'Success';
                        TestData.success         = true;
                    else
                        %=====================================================
                        % CASE: Actual and expected exception classes are NOT
                        %       identical.
                        %=====================================================
                        TestData.resultDescrText = ...
                            ['Test expected exception, and generated an', ...
                            ' actual exception, but the exception types differ.'];
                        TestData.success         = false;
                    end
                else
                    %================================
                    % CASE: Did not expect exception
                    %================================
                    TestData.resultDescrText = 'Test generated unexpected exception.';
                    TestData.success         = false;
                end
            else
                %===========================
                % CASE: No actual exception
                %===========================
                
                if isempty(obj.expExceptionName)
                    %================================
                    % CASE: Did not expect exception
                    %================================
                    if obj.Settings.equalsFunc(obj.expOutput, actOutput)
                        %========================================
                        % CASE: Actual result == expected result
                        %========================================
                        TestData.resultDescrText = 'Success';
                        TestData.success         = true;
                    else
                        %========================================
                        % CASE: Actual result != expected result
                        %========================================
                        
                        % NOTE: Only the location of the "FIRST" difference?!
                        [equals_recursive_equal, diffLoc, diffMsg] = ...
                            EJ_library.utils.equals_recursive(...
                            obj.expOutput, actOutput);
                        
                        if equals_recursive_equal
                            warning(...
                                ['EJ_library.utils.equals_recursive', ...
                                ' can not find the difference.'])
                        end
                        TestData.Result.diff.location = diffLoc;
                        TestData.Result.diff.message  = diffMsg;
                        TestData.resultDescrText      = ...
                            ['Actual function result differs from expected', ...
                            ' function result.'];
                        TestData.success              = false;
                        
                        if obj.Settings.keyboardOnUnequal
                            % To make it more convenient for the user.
                            expOutput = obj.expOutput;
                            keyboard
                        end
                    end
                else
                    %==========================
                    % CASE: Expected exception
                    %==========================
                    TestData.resultDescrText = ...
                        'Test did not generate an exception as expected.';
                    TestData.success         = false;
                end
                
            end

            % NOTE: .Parameters only useful for caller because class properties
            % are PUBLIC.
            TestData.Parameters          = obj;
            TestData.Result.actOutput    = actOutput;
            TestData.Result.actException = actException;
        end
    end
    
end
