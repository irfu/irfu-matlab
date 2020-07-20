function CURRENT_zv_to_current___ATEST

    % Submit data series as 2D vector.
    % Good for testing algorithm, but not input format.
    newTest2    = @(t1_zvIBIASx1, t2_zvIBIASx2, duplicatesAnomaly) (EJ_library.atest.CompareFuncResult(...
        @EJ_library.so.CURRENT_zv_to_current, ...
        {t1_zvIBIASx1(:,1), t1_zvIBIASx1(:,2)}, ...
        {t2_zvIBIASx2(:,1), t2_zvIBIASx2(:,2), duplicatesAnomaly}));
    
    newTestExc2 = @(t1_zvIBIASx1, t2) (EJ_library.atest.CompareFuncResult(...
        @EJ_library.so.CURRENT_zv_to_current, ...
        {t1_zvIBIASx1(:,1), t1_zvIBIASx1(:,2)}, 'MException'));
    
    
    
    tl = {};
    
    tl{end+1} = newTest2([...
        0,  10; ...
        1,  11; ...
        2,  12; ...
        ], [ ...
        0, 10; ...
        1, 11; ...
        2, 12], ...
        0);
    
    tl{end+1} = newTest2([...
        0, NaN; ...
        1,  11; ...
        2,  12; ...
        3, NaN; ...
        4   14 ...
        ], [ ...
        1, 11; ...
        2, 12; ...
        4, 14], ...
        0);
    
    % Duplicate anomaly
    tl{end+1} = newTest2([...
        0, NaN; ...
        1,  11; ...
        2,  12; ...
        2,  12; ...
        2,  12; ...
        3, NaN; ...
        4,  14; ...
        4,  14; ...
        5,  15 ...
        ], [ ...
        1, 11; ...
        2, 12; ...
        4, 14; ...
        5, 15], ...
        1);
    
    % Illegal duplicates
    tl{end+1} = newTestExc2([...
        0, NaN; ...
        1,  11; ...
        2,  12; ...
        2,  12; ...
        2,  12.1; ...
        3, NaN; ...
        4,  14; ...
        4,  14; ...
        5,  15 ...
        ]);
    
    %tl{end+1} = newTest();
    %tl{end+1} = newTest();
    %tl{end+1} = newTest();
    %tl{end+1} = newTest();
    %tl{end+1} = newTest();
    %tl{end+1} = newTest();
    
    EJ_library.atest.run_tests(tl)
end