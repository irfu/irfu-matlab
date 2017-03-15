% dm_utils_TEST - Automated test code for dm_utils.
% 
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%

% PROPOSAL: Generic function for testing a function without side effects. User submits function pointer, list of
% argument lists, list of expected return results.
%   NOTE: Has to handle approximate numeric results.

function dm_utils_TEST
    %find_last_same_sequence_TEST
    %convert_N_to_1_SPR_ACQUISITION_TIME_TEST
    %convert_N_to_1_SPR_Epoch_TEST
    %convert_N_to_1_SPR_redistribute_TEST
    %convert_N_to_1_SPR_repeat_TEST
    find_sequences_TEST
end



function find_last_same_sequence_TEST

    function i = split_into_sequences(varargin)
	% Utility function
        i = [];
        i_first = 1;        
        while i_first <= length(varargin{1})
            i(end+1) = i_first;
            i_last = bicas.dm_utils.find_last_same_sequence(i_first, varargin{:});            
            i_first = i_last + 1;
        end
        i(end+1) = i_first;
    end

    i_res = {};
    i_exp = {};
    i_res{end+1} = split_into_sequences([1,1,1]');
    i_exp{end+1} = [1,4];
    i_res{end+1} = split_into_sequences([1,5]');
    i_exp{end+1} = [1,2,3];
    i_res{end+1} = split_into_sequences([1,1,1,5,6,6]');
    i_exp{end+1} = [1,4,5,7];
    i_res{end+1} = split_into_sequences([1,1,1,5,6,6]', [2,2,2,7,4,4]');
    i_exp{end+1} = [1,4,5,7];
    i_res{end+1} = split_into_sequences([1,1,5,5,6,6]', [2,2,2,7,4,4]');
    i_exp{end+1} = [1,3,4,5,7];
    i_res{end+1} = split_into_sequences([NaN]');
    i_exp{end+1} = [1,2];
    i_res{end+1} = split_into_sequences([NaN,NaN]');
    i_exp{end+1} = [1,3];
    i_res{end+1} = split_into_sequences([NaN,1,1,5,6,6]');
    i_exp{end+1} = [1,2,4,5,7];
    i_res{end+1} = split_into_sequences([1,1,NaN,6,6]');
    i_exp{end+1} = [1,3,4,6];
    i_res{end+1} = split_into_sequences([1,1,6,6,NaN]');
    i_exp{end+1} = [1,3,5,6];
    
    for k = 1:length(i_res)
        if ~bicas.utils.equals_tolerance(i_res{k}, i_exp{k}, 0)
            i_res{k}
            i_exp{k}
            error('FAIL')
        end
    end
end



function find_sequences_TEST
    % Rewrite to test expected consistency?
    % iFirst(2:end) + 1 == iLast(1:end-1)
    % numel(unique(A(iFirst(i):iLast(i)))) == 1
    % A(iLast(1:end-1)) ~= A(iFirst(2:end))    // if only one vector.

    args = {};
    exp = {};
    args{end+1} = {zeros(0,1)};
    exp{end+1} = zeros(0,2);
    args{end+1} = {[1,1,1]'};
    exp{end+1} = [1,3];
    args{end+1} = {[1,5]'};
    exp{end+1} = [1,1; 2,2];
    args{end+1} = {[1,1,1,5,6,6]'};
    exp{end+1} = [1,3; 4,4; 5,6];
    args{end+1} = {[1,1,1,5,6,6]', [2,2,2,7,4,4]'};
    exp{end+1} = [1,3;4,4;5,6];
    args{end+1} = {[1,1,5,5,6,6]', [2,2,2,7,4,4]'};
    exp{end+1} = [1,2; 3,3; 4,4; 5,6];
    args{end+1} = {[1,1, NaN,NaN, 6,6,6]', [2,2,2, 7, 4,4,4]'};
    exp{end+1} = [1,2; 3,3; 4,4; 5,7];
    args{end+1} = {[NaN]'};
    exp{end+1} = [1,1];
    args{end+1} = {[NaN,NaN]'};
    exp{end+1} = [1,2];
    args{end+1} = {[NaN,1,1,5,6,6]'};
    exp{end+1} = [1,1; 2,3; 4,4; 5,6];
    args{end+1} = {[1,1,NaN,6,6]'};
    exp{end+1} = [1,2; 3,3; 4,5];
    args{end+1} = {[1,1,6,6,NaN]'};
    exp{end+1} = [1,2; 3,4; 5,5];
    
    for k = 1:length(args)
        [res1, res2] = bicas.dm_utils.find_sequences(args{k}{:});
        res = [res1(:), res2(:)];
        if ~bicas.utils.equals_tolerance(res, exp{k}, 0)
            args{k}{:}
            exp{k}
            res
            error('FAIL')
        end
    end
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

