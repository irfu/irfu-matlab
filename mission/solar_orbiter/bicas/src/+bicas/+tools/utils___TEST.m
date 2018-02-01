% utils_TEST - Automated test code for utils.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-03-15
%

function utils_TEST
    select_array_structs_TEST
end


function select_array_structs_TEST

    args = {};
    exp = {};
    args{end+1} = {struct('a', {'qwe'}), 'a', {NaN}};
    exp{end+1}  = {struct('a', {}), [0]};
    args{end+1} = {struct('a', {'qwe', 123, NaN}), 'a', {'qwe', 123}};
    exp{end+1}  = {struct('a', {'qwe', 123}), [1 1 0]'};
    
    args{end+1} = {struct('a', {1 1 2 2 3 3}, 'b', {5 6 5 6 5 6}), 'a', {2 4}, 'b', {5}};
    exp{end+1}  = {struct('a', {2}, 'b', {5}), [0 0 1 0 0 0]'};
    
    for k = 1:length(args)
        [res1, res2] = bicas.utils.select_array_structs(args{k}{:});
        res = {res1, res2};
        if ~isequaln(res, exp{k})    % NOTE: isequaln considers 0==false, 1==true.
            a = args{k}
            e = exp{k}
            r = res
            error('FAIL')
        end
    end
end
