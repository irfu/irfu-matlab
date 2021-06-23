function interpolate_nearest___ATEST
    newTest = @(xDist, x1, y1, x2, outputOrExc) (EJ_library.atest.CompareFuncResult...
        (@bicas.utils.interpolate_nearest, ...
        {xDist, x1, y1, x2}, outputOrExc));
    
    tl = {};
    tl{end+1} = newTest(0,   [],    [],      [1 2 3 4 5],   {[NaN NaN NaN NaN NaN]});
    
    tl{end+1} = newTest(0,   [3],   [13],    [1 2 3 4 5],   {[NaN NaN 13  NaN NaN]});
    tl{end+1} = newTest(1,   [3],   [13],    [1 2 3 4 5],   {[NaN 13  13  13  NaN]});
    tl{end+1} = newTest(2,   [3],   [13],    [1 2 3 4 5],   {[13  13  13  13  13]});
    
    tl{end+1} = newTest(0,   [3 4], [13 14], [1 2 3 4 5 6], {[NaN NaN 13 14 NaN NaN]});
    tl{end+1} = newTest(0.9, [3 4], [13 14], [1 2 3 4 5 6], {[NaN NaN 13 14 NaN NaN]});
    tl{end+1} = newTest(1,   [3 4], [13 14], [1 2 3 4 5 6], {[NaN 13  13 14 14  NaN]});
    tl{end+1} = newTest(Inf, [3 4], [13 14], [1 2 3 4 5 6], {[13  13  13 14 14  14 ]});
    
    tl{end+1} = newTest(1,   [3 4], [13 14], [], {[]});
    
    tl{end+1} = newTest(1,   [3 4], [13 14], [3 4],  {[13 14]});
    tl{end+1} = newTest(1,   [3 4], [13 14], [3 4]', {[13 14]'});
    
    tl{end+1} = newTest(1,   [], [], [], {[]});
    tl{end+1} = newTest(1,   [], [], ones(13,0), {ones(13,0)});
    
    EJ_library.atest.run_tests(tl)
end