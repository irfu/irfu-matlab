function zv_TC_to_current___ATEST
    newTest    = @(t1, zvIBIASx1, t2, zvIBIASx2, mitigatedDuplicates) (EJ_library.atest.CompareFuncResult(...
        @EJ_library.so.zv_TC_to_current, {t1, zvIBIASx1, t2}, {zvIBIASx2, mitigatedDuplicates}));
    newTestExc = @(t1, zvIBIASx1, t2) (EJ_library.atest.CompareFuncResult(@function_to_test, {t1, zvIBIASx1, t2}, 'MException'));
    
    tl = {};
    
    % General tests. Includes values before & after data time range.
    % Duplicate anomaly.
    tl{end+1} = newTest([0   1   2   3 4 5], [11    nan     nan     12 nan nan], [-1:6], [NaN 11 11 11 12 12 12 12], 0);
    % No duplicate anomaly.
    tl{end+1} = newTest([0 0 1 1 2 2 3 4 5], [11 11 nan nan nan nan 12 nan nan], [-1:6], [NaN 11 11 11 12 12 12 12], 1);
    
    % Duplicate timestamps, but not duplicate bias.
    tl{end+1} = newTestExc([0 0 1   2   3 4 5], [11 10 nan nan nan nan 12 nan nan], [-1:6]);    
    % Sorted timestamps for the antenna, but not locally.
    %tl{end+1} = newTest(   [1   0   2   3 4 5], [11    nan     nan     12 nan nan], [-1:6], [nan nan 11 11 12 12 12 12], 0);
    tl{end+1} = newTestExc([1   0   2   3 4 5], [11    nan     nan     12 nan nan], [-1:6]);
    
    %tl{end+1} = 
    %tl{end+1} = 
    %tl{end+1} = 
    %tl{end+1} = 
    %tl{end+1} = 
    
    EJ_library.atest.run_tests(tl)
end
