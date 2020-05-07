%
% Automatic test code.
%
% NOTE: Can not test very much due to how EJ_library.utils.tabulated_transform is built now (no evalutation). Test
% exists for ~historical reasons.
%
function tabulated_transform___ATEST
    % PROPOSAL: Test toward_zero_at_high_freq
    
    newTest_constructor_crash = @(input, expOutputOrExcName) (EJ_library.atest.CompareFuncResult(...
        @test_constructor_crash, ...
        input, ...
        expOutputOrExcName, ...
        'equalsFunc', @(a,b) (EJ_library.utils.equals_recursive(a,b, 'epsilon', 1e-7))));
    
    
    
    tl = {};
    
    tl{end+1} = newTest_constructor_crash({[]}, 'MException');
    tl{end+1} = newTest_constructor_crash({[], []}, 'MException');
    tl{end+1} = newTest_constructor_crash({[1,2], [1,2]}, {});
    tl{end+1} = newTest_constructor_crash({[1,2], [1]}, 'MException');

    
    
    EJ_library.atest.run_tests(tl)
end



function test_constructor_crash(varargin)
    Transf = EJ_library.utils.tabulated_transform(varargin{:});
end
