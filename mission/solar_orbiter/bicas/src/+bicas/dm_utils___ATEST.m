% dm_utils_TEST - Automated test code for functions in dm_utils. Does not necessarily test all functions.
%
% NOTE: Does NOT test all functions in dm_utils.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%

% PROPOSAL: Generic function for testing a function without side effects. User submits function pointer, list of
% argument lists, list of expected return results.
%   NOTE: Has to handle approximate numeric results.

function dm_utils___ATEST
    find_sequences_TEST
    convert_N_to_1_SPR_ACQUISITION_TIME_TEST
    convert_N_to_1_SPR_Epoch_TEST
    convert_N_to_1_SPR_redistribute_TEST
    convert_N_to_1_SPR_repeat_TEST
end



function find_sequences_TEST
    
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.dm_utils.find_sequences, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[]', []'}, 'MException');    % NOTE: size([]) = 0x0 ==> Not column vector
    tl{end+1} = new_test({ones(0,1), ones(0,1)}, {[], []});
    tl{end+1} = new_test({[1]'                          }, {[1], [1]});
    tl{end+1} = new_test({[1]',           [3]'          }, {[1], [1]});
    tl{end+1} = new_test({[1,1,1]',       [3,3,3]'      }, {[1], [3]});    
    tl{end+1} = new_test({[1,1,1,2,2,2]', [3,3,3,4,4,4]'}, {[1,4], [3,6]});
    
    EJ_library.atest.run_tests(tl)
end



function convert_N_to_1_SPR_ACQUISITION_TIME_TEST
% NOTE: Indirectly tests convert_N_to_1_SPR_Epoch too.

    % NOTE: Tests should actually be independent of the exact value(!)
    ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];

    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.dm_utils.convert_N_to_1_SPR_ACQUISITION_TIME, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test( {uint32([123, 100]), 1, 100, ACQUISITION_TIME_EPOCH_UTC}, {uint32([123, 100])});    
    tl{end+1} = new_test(...
        {uint32([123, 100; 4, 65486]),                   4, [65536/100; 65536/100], ACQUISITION_TIME_EPOCH_UTC}, ...
        {uint32([123, 100; 123, 200; 123, 300; 123, 400; 4, 65486; 5, 50; 5, 150; 5, 250])});
    
    EJ_library.atest.run_tests(tl)
end



function convert_N_to_1_SPR_Epoch_TEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.dm_utils.convert_N_to_1_SPR_Epoch, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({int64(1000),         3, 1e6},        {int64([1000; 2000; 3000])});
    tl{end+1} = new_test({int64([2000; 3000]), 1, [1e7; 1e7]}, {int64([2000; 3000])  });
    tl{end+1} = new_test({int64([2000; 3000]), 2, [1e7; 1e8]}, {int64([2000; 2100; 3000; 3010])});
    
    EJ_library.atest.run_tests(tl);
end



function convert_N_to_1_SPR_redistribute_TEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.dm_utils.convert_N_to_1_SPR_redistribute, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[1,2,3; 4,5,6; 7,8,9; 10,11,12]}, {(1:12)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]},       {(1:11)'});
    tl{end+1} = new_test({[1,2,3,4,5,6,7,8,9,10,11]'},      {(1:11)'});
    
    EJ_library.atest.run_tests(tl);
end



function convert_N_to_1_SPR_repeat_TEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.dm_utils.convert_N_to_1_SPR_repeat, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({[5;6], 2},      {[5;5;6;6]});    
    tl{end+1} = new_test({[5], 2},        {[5;5]});    
    tl{end+1} = new_test({zeros(0,1), 2}, {zeros(0,1)});
    
    EJ_library.atest.run_tests(tl);
end
