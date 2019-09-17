% Automated test code for functions in proc_utils.
%
% NOTE: Does NOT test all functions in proc_utils.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%

% PROPOSAL: Generic function for testing a function without side effects. User submits function pointer, list of
% argument lists, list of expected return results.
%   NOTE: Has to handle approximate numeric results.

function proc_utils___ATEST
    find_sequences___ATEST
    convert_N_to_1_SPR_ACQUISITION_TIME___ATEST
    convert_N_to_1_SPR_Epoch___ATEST
    convert_N_to_1_SPR_redistribute___ATEST
    convert_N_to_1_SPR_repeat___ATEST
    set_NaN_after_snapshots_end___ATEST
end



function find_sequences___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.find_sequences, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[]', []'}, 'MException');    % NOTE: size([]) = 0x0 ==> Not column vector
    tl{end+1} = new_test({ones(0,1), ones(0,1)}, {[], []});
    tl{end+1} = new_test({[1]'                          }, {[1], [1]});
    tl{end+1} = new_test({[1]',           [3]'          }, {[1], [1]});
    tl{end+1} = new_test({[1,1,1]',       [3,3,3]'      }, {[1], [3]});    
    tl{end+1} = new_test({[1,1,1,2,2,2]', [3,3,3,4,4,4]'}, {[1,4], [3,6]});
    
    EJ_library.atest.run_tests(tl)
end



function convert_N_to_1_SPR_ACQUISITION_TIME___ATEST
% NOTE: Indirectly tests convert_N_to_1_SPR_Epoch too.

    % NOTE: Tests should actually be independent of the exact value(!)
    ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];

    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.convert_N_to_1_SPR_ACQUISITION_TIME, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test( {uint32([123, 100]), 1, 100, ACQUISITION_TIME_EPOCH_UTC}, {uint32([123, 100])});    
    tl{end+1} = new_test(...
        {uint32([123, 100; 4, 65486]),                   4, [65536/100; 65536/100], ACQUISITION_TIME_EPOCH_UTC}, ...
        {uint32([123, 100; 123, 200; 123, 300; 123, 400; 4, 65486; 5, 50; 5, 150; 5, 250])});
    
    EJ_library.atest.run_tests(tl)
end



function convert_N_to_1_SPR_Epoch___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.convert_N_to_1_SPR_Epoch, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({int64(1000),         3, 1e6},        {int64([1000; 2000; 3000])});
    tl{end+1} = new_test({int64([2000; 3000]), 1, [1e7; 1e7]}, {int64([2000; 3000])  });
    tl{end+1} = new_test({int64([2000; 3000]), 2, [1e7; 1e8]}, {int64([2000; 2100; 3000; 3010])});
    
    EJ_library.atest.run_tests(tl);
end



function convert_N_to_1_SPR_redistribute___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.convert_N_to_1_SPR_redistribute, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[1,2,3; 4,5,6; 7,8,9; 10,11,12]}, {(1:12)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]},       {(1:11)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]'},      {(1:11)'});
    
    EJ_library.atest.run_tests(tl);
end



function convert_N_to_1_SPR_repeat___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.convert_N_to_1_SPR_repeat, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[5;6], 2},      {[5;5;6;6]});    
    tl{end+1} = new_test({[5], 2},        {[5;5]});    
    tl{end+1} = new_test({zeros(0,1), 2}, {zeros(0,1)});
    
    EJ_library.atest.run_tests(tl);
end



function set_NaN_after_snapshots_end___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.set_NaN_after_snapshots_end, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({ones(0,4),              ones(0,1)},   {ones(0,4)});
    tl{end+1} = new_test({[0,1,2],                [3]},   {[0,1,2]});
    tl{end+1} = new_test({[0,1,2,3,4; 5,6,7,8,9], [2;4]}, {[0,1,NaN,NaN,NaN; 5,6,7,8,NaN]});
    
    EJ_library.atest.run_tests(tl);
end
