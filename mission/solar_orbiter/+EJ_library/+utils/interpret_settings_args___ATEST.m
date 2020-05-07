%
% Automatic test code for EJ_library.utils.interpret_settings_args.
%
function interpret_settings_args___ATEST
    
    % NOTE: Does not check recursive behaviour.

    newTest = @(input, expOutputOrExcName) (EJ_library.atest.CompareFuncResult(...
        @EJ_library.utils.interpret_settings_args, input, expOutputOrExcName));

    tl = {};

    
    tl{end+1} = newTest({struct(),                          {}},                                       {struct()});
    tl{end+1} = newTest({struct('a', 0),                    {}},                                       {struct('a', 0)});
    tl{end+1} = newTest({struct('a', 0),                    {'a', 1}},                                 {struct('a', 1)});
    tl{end+1} = newTest({struct('enabled', 0),              {}},                                       {struct('enabled', 0)});
    tl{end+1} = newTest({struct('enabled', 0),              {'enabled', 1}},                           {struct('enabled', 1)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2),            {'b', 222}},                               {struct('a', 1, 'b', 222)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2),            {struct('b', 22)}},                        {struct('a', 1, 'b',  22)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2, 'c', 3),    {struct('b', 22), 'c', 333}},              {struct('a', 1, 'b',  22, 'c', 333)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2, 'c', 3),    {struct('b', 22), 'b', 222, 'c', 333}},    {struct('a', 1, 'b', 222, 'c', 333)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2, 'c', 3),    {'b', 222, 'c', 333}},                     {struct('a', 1, 'b', 222, 'c', 333)});
    tl{end+1} = newTest({struct('a', 1, 'b', 2, 'c', 3),    {struct(), 'b', 222, 'c', 333}},           {struct('a', 1, 'b', 222, 'c', 333)});

    % Add keys not present in DefaultSettings (permitted!).
    tl{end+1} = newTest({struct('a', 1, 'b', 2),            {'c', 333}},                               {struct('a', 1, 'b', 2,  'c', 333)});
    tl{end+1} = newTest({struct('a', 1),                    {struct('b', 22), 'c', 333}},              {struct('a', 1, 'b', 22, 'c', 333)});

    EJ_library.atest.run_tests(tl)
end
