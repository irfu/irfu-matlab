% dm_utils_TEST - Automated test code for functions in dm_utils. Does not necessarily test all functions.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%

% PROPOSAL: Generic function for testing a function without side effects. User submits function pointer, list of
% argument lists, list of expected return results.
%   NOTE: Has to handle approximate numeric results.

function dm_utils_TEST
    find_last_same_sequence_TEST
    convert_N_to_1_SPR_ACQUISITION_TIME_TEST
    convert_N_to_1_SPR_Epoch_TEST
    convert_N_to_1_SPR_redistribute_TEST
    convert_N_to_1_SPR_repeat_TEST
    find_sequences_TEST
end



function find_sequences_TEST
    % Rewrite to test expected consistency?
    % iFirst(2:end) + 1 == iLast(1:end-1)
    % numel(unique(A(iFirst(i):iLast(i)))) == 1
    % A(iLast(1:end-1)) ~= A(iFirst(2:end))    // if only one vector.



end



function convert_N_to_1_SPR_ACQUISITION_TIME_TEST
% NOTE: Indirectly tests convert_N_to_1_SPR_Epoch too.

    AT_res = {};
    AT_exp = {};
    AT_res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_ACQUISITION_TIME(uint32([123, 100]), 1, 100);
    AT_exp{end+1} = uint32([123, 100]);
    AT_res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_ACQUISITION_TIME(uint32([123, 100; 4, 65486]), 4, [65536/100; 65536/100]);
    AT_exp{end+1} = uint32([123, 100; 123, 200; 123, 300; 123, 400; 4, 65486; 5, 50; 5, 150; 5, 250]);
    
    i = 1;
    while i<=length(AT_res)
        if ~bicas.utils.equals_tolerance(AT_res{i}, AT_exp{i}, 0)
            AT_res{i}
            AT_exp{i}
            error('FAIL')
        end
        i = i + 1;
    end
    
end



function convert_N_to_1_SPR_Epoch_TEST
    args = {};
    exp = {};
    res = {};
    
    args{end+1} = {int64(1000), 3, 1e6};
    exp{end+1} = int64([1000; 2000; 3000]);
    
    args{end+1} = {int64([2000; 3000]), 1, [1e7;1e7]};
    exp{end+1} = int64([2000; 3000]);
    
    args{end+1} = {int64([2000; 3000]), 2, [1e7; 1e8]};
    exp{end+1} = int64([2000; 2100; 3000; 3010]);
    
    for k = 1:length(args)
        res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_Epoch(args{k}{:});
        
        if ~bicas.utils.equals_tolerance(res{k}, exp{k}, 0)
            res{k}
            exp{k}
            error('FAIL')
        end
    end
end



function convert_N_to_1_SPR_redistribute_TEST
    res = {};
    exp = {};
    res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_redistribute([1,2,3; 4,5,6; 7,8,9; 10,11,12]);
    exp{end+1} = (1:12)';
    res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_redistribute([1,2,3,4,5,6,7,8,9,10,11]);
    exp{end+1} = (1:11)';
    res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_redistribute([1,2,3,4,5,6,7,8,9,10,11]');
    exp{end+1} = (1:11)';
    
    for k = 1:length(res)
        if ~bicas.utils.equals_tolerance(res{k}, exp{k}, 0)
            res{k}
            exp{k}
            error('FAIL')
        end
    end
end



function convert_N_to_1_SPR_repeat_TEST
    args = {};
    exp = {};
    res = {};
    args{end+1} = {[5;6], 2};
    exp{end+1} = [5;5;6;6];
    args{end+1} = {[5], 2};
    exp{end+1} = [5;5];
    args{end+1} = {zeros(0,1), 2};
    exp{end+1} = zeros(0,1);
    
    for k = 1:length(args)
        res{end+1} = bicas.dm_utils.convert_N_to_1_SPR_repeat(args{k}{:});
        
        if ~bicas.utils.equals_tolerance(res{k}, exp{k}, 0)
            res{k}
            exp{k}
            error('FAIL')
        end
    end
end

