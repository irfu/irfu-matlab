%
% Function for running automatic tests of code. Accepts list of objects, calls a
% test method in every object, and then displays return result from every such
% call.
%
% Stops after first failure.
%
%
% FUNCTION HANDLES FOR TESTING CLASS METHODS
% ==========================================
% Class instance (not static?) methods (not constructors):
%   @methodNameExcludingClassOrPackages, and having the first input argument
%   interpreted as the object/class instance.
% Class instance constructor
%   @classNameIncludingPackages
%
%
% ARGUMENTS
% =========
% testList : Cell array of objects.
%       Every objects must have a method TestData = run().
%           TestData : Struct with fields
%               .success     : Boolean
%               .resultDescr
%                   Plain English text string describing how the test failed
%                   (there can be multiple ways).
%               .Parameters  : Arbitrary, but likely struct.
%                   The input values, which define the test.
%               .Result      : Arbitrary, but likely struct.
%                   The output values, defined by the test result.
%
%
% Initially created 2019-01-21 by Erik P G Johansson.
%
function run_tests(testList)
    % PROPOSAL: Setting quitOnFirstFail.
    % PROPOSAL: Setting triggerKeyboard.
    % PROPOSAL: Setting ~printResultDepth.
    %
    % PROPOSAL: Optional whether to stop on first failure.
    %
    % PROPOSAL: On failure, print function call that generated failure so that user can call manually from CLI.
    %   NOTE: Can not alway generate usable call.
    %       Ex: Can not generate package path.
    %       Ex: Can not call inner or nested function from CLI.
    
    Settings.printParentSeparately = true;
    Settings.stringsEscape         = true;
    Settings.stringsSsMaxLength    = 120;
%    Settings.printParentSeparately = false;   % Does not work due to bug?

    assert(iscell(testList), 'Argument testList is not a cell array.')
    
    for iTest = 1:numel(testList)
        Test = testList{iTest};
        assert(isobject(Test) & isscalar(Test), ...
            '"Test" is not a scalar object. Argument "testList" is misconfigured.')
        
        TestData = Test.run();
        EJ_library.assert.struct(TestData, ...
            {'success', 'resultDescrText', 'Parameters', 'Result'}, {});
        
        fprintf('TEST %2i: ', iTest)
        if TestData.success
            fprintf('OK\n')
        else
            fprintf('FAILED: %s\n', TestData.resultDescrText)
            EJ_library.utils.print_variable_recursively('    Parameters', TestData.Parameters, Settings);
            EJ_library.utils.print_variable_recursively('    Result',     TestData.Result,     Settings);
            %keyboard
            
            break
        end
    end
end
