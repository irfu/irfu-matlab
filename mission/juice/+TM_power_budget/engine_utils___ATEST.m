function engine_utils___ATEST
    % function [xp, y1, yArray] = linear_extrapolate_limit(x0, x1, y0, dydx, xArray, yArray, yMin, yMax)
    args = {};
    exp = {};
    %args{end+1} = {0, 1, 0, 1, [], [], -1, 2};
    %exp{end+1}  = {1, 1, zeros(0,1)};
    %args{end+1} = {10, 12, 0, 1, [], [], 0, 1};
    %exp{end+1}  = {11, 2, zeros(0,1)};
    %args{end+1} = {10, 12, 1, 1, [], [], 0, 1};
    %exp{end+1}  = {10, 3, zeros(0,1)};
    %args{end+1} = {10, 10, 1, 1, [], [], 0, 1};
    %exp{end+1}  = {10, 1, zeros(0,1)};
    args{end+1} = {22933.3333333333, 43200, 3.5527136788005e-14, -0.025, [], [], 0, Inf};   % May produce xp==x0, yp~=y0
    exp{end+1}  = {22933.3333333333, 0, zeros(0,1)};
    %args{end+1} = {};
    %exp{end+1}  = {};
    
    for k = 1:length(args)
        [res1, res2, res3] = TM_power_budget.engine_utils.linear_extrapolate_limit(args{k}{:});
        res = {res1(:), res2(:), res3(:)};
        if ~isequal(res, exp{k})
            args{k}{:}
            exp{k}
            res
            error('FAIL')
        end
    end
end