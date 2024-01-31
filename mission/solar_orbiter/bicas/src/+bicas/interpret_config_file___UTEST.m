%
% matlab.unittest automatic test code for bicas.interpret_config_file().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-10, from older test code first created 2018-01-25.
%
classdef interpret_config_file___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Add tests for illegal config files.



  properties(Constant)
    L = bicas.Logger('none', false);
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      % One output variable.
      function test(configFileRowsCa, expOutput)
        actOutput = bicas.interpret_config_file(configFileRowsCa, testCase.L);
        testCase.verifyEqual(actOutput, expOutput)
      end
      %===================================================================

      test({'# Comment'},               containers.Map('KeyType', 'char', 'ValueType', 'char'));
      test({'key="value"'},             containers.Map({'key'}, {'value'}));
      test({'key="value"   # Comment'}, containers.Map({'key'}, {'value'}));
      test({...
        '# Comment', ...
        '', ...
        '   ', ...
        'key.1="value1"', ...
        'key_2="value2"   # Comment', ...
        'key-3  =   ""' ...
        }, containers.Map({'key.1', 'key_2', 'key-3'}, {'value1', 'value2', ''}));

      test({...
        'key_1 = "value1"', ...
        'key_2 = "value2"   # Comment', ...
        'key_1 = "value1_new"' ...
        }, containers.Map({'key_1', 'key_2'}, {'value1_new', 'value2'}));

    end



  end    % methods(Test)



end
