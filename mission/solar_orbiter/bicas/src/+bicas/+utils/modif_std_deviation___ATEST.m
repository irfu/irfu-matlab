%
% Automatic test code
%
function modif_std_deviation___ATEST
    
    tl = {};
    
    %function mstd = modif_std_deviation(v, ref, iDim)
    
    % Empty data.
    tl = add_test(tl, zeros(0,3), 5, 1, NaN(1,3));
    tl = add_test(tl, zeros(0,3), 5, 2, NaN(0,1));
    tl = add_test(tl, zeros(0,3), 5, 3, NaN(0,3,1));
    
    % 1D vector.
    tl = add_test(tl, ones(1,3), 5, 1, NaN( 1,3,1));
    tl = add_test(tl, ones(1,3), 5, 2, ones(1,1,1) * sqrt(3*4^2 / 2));
    tl = add_test(tl, ones(1,3), 5, 3, NaN( 1,3,1));
    
    % One data point.
    tl = add_test(tl, [0], 0, 1, NaN);
    tl = add_test(tl, [5], 5, 1, NaN);
    tl = add_test(tl, [5], 4, 1, NaN);
    tl = add_test(tl, [5], 4, 2, NaN);
    
    tl = add_test(tl, [5,6,7; 2,3,4], 4, 1, sqrt([(1^2+2^2),     2^2+1^2,    3^2+0^2] / 1));
    tl = add_test(tl, [5,6,7; 2,3,4], 4, 2, sqrt([(1^2+2^2+3^2); 2^2+1^2+0^2]         / 2));
    tl = add_test(tl, [5,6,7; 2,3,4], 4, 3, NaN(2,3,1));
    
    % Test NaN sample.
    tl = add_test(tl, [NaN,6,7; 2,3,4], 4, 1, sqrt([NaN, 2^2+1^2,     3^2+0^2] / 1));
    tl = add_test(tl, [NaN,6,7; 2,3,4], 4, 2, sqrt([NaN; 2^2+1^2+0^2]          / 2));
    tl = add_test(tl, [NaN,6,7; 2,3,4], 4, 3, NaN(2,3,1));
    
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    %tl = add_test(tl, );
    
    EJ_library.atest.run_tests(tl)
end



function tl = add_test(tl, v, ref, iDim, mstd)
    
    tl{end+1} = EJ_library.atest.CompareFuncResult(...
        @bicas.utils.modif_std_deviation, ...
        {v, ref, iDim}, ...
        {mstd});
end



% function tl = add_test_exc(tl, inputs)
%     assert(iscell(inputs))
%     
%     tl{end+1} = EJ_library.atest.CompareFuncResult(...
%         @bicas.utils.modif_std_deviation, ...
%         inputs, ...
%         'MException');
% end
