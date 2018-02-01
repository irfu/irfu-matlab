% equals_tolerance_TEST - Automated test code for equals_tolerance.
% 
% NOTE: Currently does not test for epsilon<>0.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
function equals_tolerance_TEST

% PROPOSAL: Automatically test symmetry of A, B.

    arguments = {};
    result_expected = {};
    
    
    
    % Construct list of unique test values. Unique ==> List index can be used to determine equality.
    % NOTE: Only strictly works for epsilon == 0 but one can of course choose values so that they are unique for a
    % non-zero epsilon too.
    A = [1,2,3; 4,   5,6];
    B = [1,2,3; 4,-Inf,6];
    C = [1,2,3; 4, Inf,6];
    D = [1,2,3; 4, NaN,6];    
    
    E1(:,:,1) = [ 1,2; 3,4];
    E1(:,:,2) = [ 1,2; 3,4] + 10;
    E2 = E1;
    E2(1,1,1) = 12;
    
    testValues = {};
    testValues = [testValues, E1, E2];
    testValues = [testValues, A, B, C, D];
    testValues = [testValues, zeros(0,0), zeros(0,1), zeros(1,0)];
    testValues = [testValues, zeros(2,3), zeros(3,2)];
    testValues = [testValues, reshape(11:16, 2,3), reshape(11:16, 3,2)];
    
    % Create tests for all combinations of test values.    
    for i = 1:length(testValues)
        for j = 1:length(testValues)    
            % NOTE: For loops ==> Will check i vs j AND j vs i.
            arguments      {end+1} = {testValues{i}, testValues{j}, 0};
            result_expected{end+1} = (i==j);
        end
    end
    
    

    for iTest = 1:length(arguments)
        result_actual = bicas.utils.equals_tolerance(arguments{iTest}{:});
        
        if result_actual == result_expected{iTest}    % NOTE: Do not want to use equals_tolerance for this test.
            disp('TEST OK')
        else
            % Print values
            A = arguments{iTest}{1}
            B = arguments{iTest}{2}
            result_actual
            result_exp = result_expected{iTest}
            
            error('TEST FAILED')
        end
    end
end
