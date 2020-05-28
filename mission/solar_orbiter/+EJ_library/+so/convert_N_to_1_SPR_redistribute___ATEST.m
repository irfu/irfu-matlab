function convert_N_to_1_SPR_redistribute___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@EJ_library.so.convert_N_to_1_SPR_redistribute, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[1,2,3; 4,5,6; 7,8,9; 10,11,12]}, {(1:12)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]},       {(1:11)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]'},      {(1:11)'});
    
    A = permute(reshape([1:24], 4,3,2), [2,1,3]);   % 3D matrix.
    tl{end+1} = new_test({A},   {[(1:12)', (13:24)']});
    
    EJ_library.atest.run_tests(tl);
end
