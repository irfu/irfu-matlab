function convert_N_to_1_SPR_Epoch___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@EJ_library.so.convert_N_to_1_SPR_Epoch, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({int64(1000),         3, 1e6},        {int64([1000; 2000; 3000])});
    tl{end+1} = new_test({int64([2000; 3000]), 1, [1e7; 1e7]}, {int64([2000; 3000])  });
    tl{end+1} = new_test({int64([2000; 3000]), 2, [1e7; 1e8]}, {int64([2000; 2100; 3000; 3010])});
    
    EJ_library.atest.run_tests(tl);
end
