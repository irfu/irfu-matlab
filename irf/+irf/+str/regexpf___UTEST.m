%
% matlab.unittest automatic test code for irf.str.regexpf().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16, converted from other test code.
%
classdef regexpf___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      function test(inputsCa, expOutputsCa)
        actOutputs = cell(size(expOutputsCa));
        [actOutputs{:}] = irf.str.regexpf(inputsCa{:});
        testCase.verifyEqual(actOutputs, expOutputsCa)
      end
      function test_exc(varargin)
        testCase.verifyError(...
          @() irf.str.regexpf(varargin{:}), ...
          ?MException)
      end
      %===================================================================

      % Empty strings (i.e. 0x0).
      % NOTE: Not sure what is desired behaviour. Just verifies de facto
      % behaviour.
      test({'', ''}, {true})

      % Illegal inputs ([] ~"empty string")
      test_exc([], 'asd')
      test_exc('asd', [])

      test({'asd',   'asd' }, {true });
      test({'asdF',  'asd' }, {false});
      test({'asd',   'asdF'}, {false});
      test({'asd',  '^asd$'}, {true });
      test({'asdF', '^asd$'}, {false});

      % Test 2D matrix IN --> 2D matrix OUT.
      test({'asd', {'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'}},          {logical([1,1,0; 0,0,0])});
      test({{'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'}, 'qwe'},          {logical([0,0,0; 1,1,0])});
      test({{'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'}, '[Qq][Ww][Ee]'}, {logical([0,0,0; 1,1,1])});

      % Scalar + empty cell array
      test({'asd', cell(0,1)}, {false(0,1)})
      test({cell(1,0), 'asd'}, {false(1,0)})

      % Test non-scalar + non-scalar IN --> Exception
      test_exc({'asd', 'qwe'}, {'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'})
      test_exc(cell(0,1), {'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'})
      test_exc({'asd', 'asd', 'ASD'; 'qwe', 'qwe', 'QWE'}, cell(0,1))
    end



  end    % methods(Test)



end
