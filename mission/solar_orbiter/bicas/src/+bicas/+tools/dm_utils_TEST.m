% dm_utils_TST - Automated test code for dm_utils.
% 
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-17
%
function dm_utils_TEST
    ACQUISITION_TIME___expand_to_sequences_TEST
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
