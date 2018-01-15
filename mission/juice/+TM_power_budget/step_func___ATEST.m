%
% Automatic test code for step_func.
%
function step_func___ATEST
% One-liner manual "test code".
% X=linspace(-5,5,1e3); [Y, YI] = TM_power_budget.step_func([-4,-2,2], [0.1,-1,2], X); plot(X,[Y; YI]); axis([min(X), max(X), min([Y,YI]), max([Y,YI])+3]); grid on

exp =  {};
args = {};

args{end+1} = {[], ...
               [], [0,   3,  8]};
exp{end+1}  = {    [NaN, NaN, NaN], ...
                   [NaN, NaN, NaN]};

args{end+1} = {[3], ...
               [2], [0,   3,  8]};
exp{end+1}  = {     [NaN, 2,  2], ...
                    [NaN, 0, 10]};
args{end+1} = {[2,    4], ...
               [1.5, -2], [0,   2,   3,    4,  9]};
exp{end+1}  = {           [NaN, 1.5, 1.5, -2, -2], ...
                          [NaN, 0,   1.5,  3, -7]};

args{end+1} = {[2,    4, 7], ...
               [1.5, -2, 1], [0,   2,   3,    4,  5,  6,  7,  9]};
exp{end+1}  = {              [NaN, 1.5, 1.5, -2, -2, -2,  1,  1], ...
                             [NaN, 0,   1.5,  3,  1, -1, -3, -1]};
                         
args{end+1} = {[2,    4, 7], ...
               [1.5, -2, 1], [3]};
exp{end+1}  = {              [1.5], ...
                             [1.5]};
                         
args{end+1} = {[2,    4, 7], ...
               [1.5, -2, 1], [0]};
exp{end+1}  = {              [NaN], ...
                             [NaN]};
                         
for k = 1:length(args)
    [result1, result2] = TM_power_budget.step_func(args{k}{:});
    if ~isequaln( exp{k}, {result1, result2} )
        k
        args{k}{1}
        args{k}{2}
        args{k}{3}
        [exp{k}{1}', exp{k}{2}']
        [result1',  result2']
        error('FAIL')
    end
    disp('TEST OK')
end

end
