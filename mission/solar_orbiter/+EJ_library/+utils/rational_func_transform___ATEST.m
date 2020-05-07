%
% Automatic test code
%
function rational_func_transform___ATEST

newTest_eval = @(input, expOutputOrExcName) (EJ_library.atest.CompareFuncResult(...
    @tf_eval, ...
    input, ...
    expOutputOrExcName));

newTest_hftz = @(input, expOutputOrExcName) (EJ_library.atest.CompareFuncResult(...
    @tf_high_freq_toward_zero, ...
    input, ...
    expOutputOrExcName));



tl = {};



tl{end+1} = newTest_eval({[],    [],    []}, 'MException');
tl{end+1} = newTest_eval({[0 0], [0 0], []}, 'MException');
tl{end+1} = newTest_eval({[1]    [2],   [0]}, {0.5});

% TEST: Non-trivial evaluation of Z(omega).
omega = [0, 10, 20; 1,2,3];
h = @(om) ((1 + 3*(i*om).^2) ./ (2 + 4*(i*om).^3));
tl{end+1} = newTest_eval({[1 0 3], [2 0 0 4], omega}, {h(omega)});



tl{end+1} = newTest_hftz({[1],     [2]  },   {false});
tl{end+1} = newTest_hftz({[1 0],   [2 1]},   {true});
tl{end+1} = newTest_hftz({[1 3],   [2 0]},   {false});
tl{end+1} = newTest_hftz({[1 3 4], [2 4 0]}, {false});
tl{end+1} = newTest_hftz({[1 3 4], [2 4 0 5]}, {true});

EJ_library.atest.run_tests(tl)

end



function Z = tf_eval(nc, dc, omega)
    Transf = EJ_library.utils.rational_func_transform(nc, dc);
    Z      = Transf.eval(omega);
end



function result = tf_high_freq_toward_zero(nc, dc)
    Transf = EJ_library.utils.rational_func_transform(nc, dc);
    result = Transf.zero_in_high_freq_limit();
end
