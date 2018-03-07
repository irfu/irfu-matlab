%
% Created 2018-02-28 by Erik Johansson, IRF Uppsala.
%
function repeat_struct_array___ATEST()
    args = {};
    exp = {};
    
    args{end+1} = {struct('t', {0,1,3}',           'x', {4,3,'qwe'}'), 't', 10, 2};   % Column vector
    exp{end+1}  =  struct('t', {0,1,3, 10,11,13}', 'x', {4,3,'qwe', 4,3,'qwe'}');
    args{end+1} = {struct('t', {0,1,3},            'x', {4,3,'qwe'}), 't', 10, 2};    % Row vector
    exp{end+1}  =  struct('t', {0,1,3, 10,11,13}', 'x', {4,3,'qwe', 4,3,'qwe'}');
    
    %args{end+1} = {struct('t', {0,1,3},            'x', {4,3,'qwe'}), 't', 10, 0};
    %exp{end+1}  =  struct('t', {{}}', 'x', {{}}');
    
    %args{end+1} = {};
    %exp{end+1}  = {};
    
    for k = 1:length(args)
        res = TM_power_budget.repeat_struct_array(args{k}{:});
        if ~isequal(res, exp{k})
            args{k}{:}
            exp{k}
            res
            error('FAIL')
        end
    end


end
