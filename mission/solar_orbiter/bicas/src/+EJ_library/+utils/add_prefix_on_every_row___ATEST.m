%
% Automatic test code for add_prefix_on_every_row.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-26
%
function add_prefix_on_every_row___ATEST

tl = {};

tl{end+1} = new_test(    'asd',      'prefix ', 'prefix asd\n');
tl{end+1} = new_test_EXC('asd\nqwe', 'prefix ');
tl{end+1} = new_test(     '',        'prefix ', 'prefix \n');

tl{end+1} = new_test('asd\n',         'prefix ', 'prefix asd\n');
tl{end+1} = new_test('asd\nqwerty\n', 'prefix ', 'prefix asd\nprefix qwerty\n');


EJ_library.atest.run_tests(tl)

end



function Test = new_test(inputStr, prefix, outputStr)

input     = {sprintf(inputStr), prefix};
expOutput = {sprintf(outputStr)};

Test = EJ_library.atest.CompareFuncResult(@EJ_library.utils.add_prefix_on_every_row, input, expOutput);

end



function Test = new_test_EXC(inputStr, prefix)

input     = {sprintf(inputStr), prefix};
expOutput = 'MException';

Test = EJ_library.atest.CompareFuncResult(@EJ_library.utils.add_prefix_on_every_row, input, expOutput);

end
