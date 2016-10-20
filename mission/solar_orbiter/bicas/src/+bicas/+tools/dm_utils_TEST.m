% dm_utils_TST - Automated test code for dm_utils.
% 
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%
function dm_utils_TEST
    ACQUISITION_TIME___expand_to_sequences_TEST
    find_last_same_sequence_TEST
    reshape_to_1_sample_per_record_TEST
end



function ACQUISITION_TIME___expand_to_sequences_TEST
% NOTE: Indirectly tests tt2000___expand_to_sequences too.

    AT_res = {};
    AT_exp = {};
    AT_res{end+1} = bicas.dm_utils.ACQUISITION_TIME___expand_to_sequences(uint32([123, 100]), 1, 100);
    AT_exp{end+1} = uint32([123, 100]);
    AT_res{end+1} = bicas.dm_utils.ACQUISITION_TIME___expand_to_sequences(uint32([123, 100; 4, 65486]), 4, [65536/100; 65536/100]);
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



function find_last_same_sequence_TEST

    function i = split_into_sequences(varargin)
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



function reshape_to_1_sample_per_record_TEST
    res = {};
    exp = {};
    res{end+1} = bicas.dm_utils.reshape_to_1_sample_per_record([1,2,3; 4,5,6; 7,8,9; 10,11,12]);
    exp{end+1} = (1:12)';
    res{end+1} = bicas.dm_utils.reshape_to_1_sample_per_record([1,2,3,4,5,6,7,8,9,10,11]);
    exp{end+1} = (1:11)';
    res{end+1} = bicas.dm_utils.reshape_to_1_sample_per_record([1,2,3,4,5,6,7,8,9,10,11]');
    exp{end+1} = (1:11)';
    
    for k = 1:length(res)
        if ~bicas.utils.equals_tolerance(res{k}, exp{k}, 0)
            res{k}
            exp{k}
            error('FAIL')
        end
    end
end

