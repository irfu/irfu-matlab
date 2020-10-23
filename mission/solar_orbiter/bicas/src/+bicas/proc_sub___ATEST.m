%
% Automatic test code
%
function proc_sub___ATEST()
    downsample_bin_sci_values___ATEST()
end



function downsample_bin_sci_values___ATEST()
    newTest    = @(zVarSegment, med, mstd) (...
        erikpgjohansson.atest.CompareFuncResult(...
        @bicas.proc_sub.downsample_bin_sci_values, ...
        {zVarSegment}, ...
        {med, mstd}));
    
    newTestExc = @(zVarSegment, med, mstd) (...
        erikpgjohansson.atest.CompareFuncResult(...
        @bicas.proc_sub.downsample_bin_sci_values, ...
        {zVarSegment}, ...
        {med, mstd}));

    ERA = zeros(1,0);
    
    tl = {};
    
    % Empty data
    tl{end+1} = newTest(zeros(0,0), ERA, ERA);
    tl{end+1} = newTest(zeros(0,2), [NaN NaN], [NaN, NaN]);
    
    % mstd=0
    tl{end+1} = newTest([1;     2;     3    ]', [1,2,3], [0,0,0]);
    tl{end+1} = newTest([1,1;   2,2;   3,3  ]', [1,2,3], [0,0,0]);
    tl{end+1} = newTest([1,1,1; 2,2,2; 3,3,3]', [1,2,3], [0,0,0]);
    
    % Average of two values (special case)
    tl{end+1} = newTest([1,2; 2,3; 3,4]',    [1.5, 2.5, 3.5], sqrt(0.5)*[1,1,1]);
    % Nominal median
    tl{end+1} = newTest([1,2,10]',    [2], sqrt( (1^2+0^2+8^2)/2 ));
    
    erikpgjohansson.atest.run_tests(tl)
end
