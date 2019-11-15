%
% Automatic test code.
%
function tabulated_transform___ATEST
% PROPOSAL: Test toward_zero_at_high_freq

newTest_eval = @(input, expOutputOrExcName) (EJ_library.atest.CompareFuncResult(...
    @tf_eval, ...
    input, ...
    expOutputOrExcName, ...
    'equalsFunc', @(a,b) (EJ_library.utils.equals_recursive(a,b, 'epsilon', 1e-7))));



tl = {};

tl{end+1} = newTest_eval({[],     [],    []}, 'MException');
tl{end+1} = newTest_eval({{[0,10], [0,1]            }, [0, 5, 10]},   {[0,0.5,1]});
tl{end+1} = newTest_eval({{[0,10], [2,3], [0,pi/2]}, [0, 5, 10]}, {[2, mean([2,3*i]), 3*i]});   % Test interpolation.
tl{end+1} = newTest_eval({{[1,10], [2,3], [0,pi/2]}, [1, 100]}, 'MException');
tl{end+1} = newTest_eval({{[1,10], [2,3], [0,pi/2]}, [0, 10]},  'MException');
tl{end+1} = newTest_eval({{[2,10], [2,3i], 'extrapolatePositiveFreqZtoZero', 1}, [0, 6, 10]}, {[2, mean([2,3i]), 3i]});   % Test interpolation.

%tl{end+1} = newTest_eval();
%tl{end+1} = newTest_eval();
%tl{end+1} = newTest_eval();

EJ_library.atest.run_tests(tl)

end



function Z = tf_eval(constructorArgs, omega)
    Transf = EJ_library.utils.tabulated_transform(constructorArgs{:});
    Z      = Transf.eval_linear(omega);
end
