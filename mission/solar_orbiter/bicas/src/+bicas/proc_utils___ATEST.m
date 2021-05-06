% Automated test code for functions in proc_utils.
%
% NOTE: Does NOT test all functions in proc_utils.
% 
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-10-17
%
function proc_utils___ATEST
% PROPOSAL: Generic function for testing a function without side effects. User submits function pointer, list of
% argument lists, list of expected return results.
%   NOTE: Has to handle approximate numeric results.
    
    set_struct_field_rows___ATEST
    convert_matrix_to_cell_array_of_vectors___ATEST
    convert_cell_array_of_vectors_to_matrix___ATEST

    % Function that is tested is commented out.
    %convert_N_to_1_SPR_ACQUISITION_TIME___ATEST
    set_NaN_after_snapshots_end___ATEST
    
end



function set_struct_field_rows___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(...
        @bicas.proc_utils.set_struct_field_rows, inputs, outputs));
    tl = {};
    
    % VERY INCOMPLETE TEST SUITE.
    
    tl{end+1} = new_test({...
        struct('asd', reshape([1:24], [4,3,2])), ...
        struct('asd', reshape([1:12], [2,3,2])), 2:3}, {...
        struct('asd', reshape([1,1,2,4, 5,3,4,8, 9,5,6,12,   13,7,8,16, 17,9,10,20, 21,11,12,24], 4,3,2)) ...
        });
    
    tl{end+1} = new_test({...
        struct('asd', [1;2;3;4;5]), ...
        struct('asd', [8;9]), [4,3]}, {...
        struct('asd', [1;2;9;8;5])});
    
%     tl{end+1} = new_test({}, {});
%     tl{end+1} = new_test({}, {});
%     tl{end+1} = new_test({}, {});
%     tl{end+1} = new_test({}, {});
    
    EJ_library.atest.run_tests(tl)
end



function convert_matrix_to_cell_array_of_vectors___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(...
        @bicas.proc_utils.convert_matrix_to_cell_array_of_vectors, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({zeros(0,1), zeros(0,1)}, {cell(0,1)});
    tl{end+1} = new_test({[1,2,3,4,5], [3]}, {{[1,2,3]}});
    tl{end+1} = new_test({[1,2,3,4,5; 6,7,8,9,0], [3; 2]}, {{[1,2,3]; [6,7]}});
    
    EJ_library.atest.run_tests(tl)
end



function convert_cell_array_of_vectors_to_matrix___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(...
        @bicas.proc_utils.convert_cell_array_of_vectors_to_matrix, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({cell(0,1),        5},   {ones(0,5), ones(0,1)});
    tl{end+1} = new_test({{[1,2,3]       }, 5},   {[1,2,3,NaN,NaN                 ], [3   ]'});
    tl{end+1} = new_test({{[1,2,3]; [1,2]}, 5},   {[1,2,3,NaN,NaN; 1,2,NaN,NaN,NaN], [3, 2]'});
    
    EJ_library.atest.run_tests(tl)
end



% function convert_N_to_1_SPR_ACQUISITION_TIME___ATEST
% % NOTE: Indirectly tests convert_N_to_1_SPR_Epoch too.
% 
%     % NOTE: Tests should actually be independent of the exact value(!)
%     ACQUISITION_TIME_EPOCH_UTC = [2000,01,01, 12,00,00, 000,000,000];
% 
%     new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.convert_N_to_1_SPR_ACQUISITION_TIME, inputs, outputs));
%     tl = {};
%     
%     tl{end+1} = new_test( {uint32([123, 100]), 1, 100, ACQUISITION_TIME_EPOCH_UTC}, {uint32([123, 100])});    
%     tl{end+1} = new_test(...
%         {uint32([123, 100; 4, 65486]),                   4, [65536/100; 65536/100], ACQUISITION_TIME_EPOCH_UTC}, ...
%         {uint32([123, 100; 123, 200; 123, 300; 123, 400; 4, 65486; 5, 50; 5, 150; 5, 250])});
%     
%     EJ_library.atest.run_tests(tl)
% end



function set_NaN_after_snapshots_end___ATEST
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_utils.set_NaN_after_snapshots_end, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({ones(0,4),              ones(0,1)},   {ones(0,4)});
    tl{end+1} = new_test({[0,1,2],                [3]},   {[0,1,2]});
    tl{end+1} = new_test({[0,1,2,3,4; 5,6,7,8,9], [2;4]}, {[0,1,NaN,NaN,NaN; 5,6,7,8,NaN]});
    
    EJ_library.atest.run_tests(tl);
end
