%
% matlab.unittest automatic test code for bicas.settings.KeyValue class.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef KeyValue___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor(T)
      Skv = bicas.settings.KeyValue(99, 'default');

      T.assertEqual(Skv.valuesCa,       {99})
      T.assertEqual(Skv.valueSourcesCa, {'default'})
    end



    function test_override(T)
      Skv = bicas.settings.KeyValue(99,  'default');
      Skv = Skv.override(          123, 'override');

      % NOTE: Column arrays.
      T.assertEqual(Skv.valuesCa,       {99;        123})
      T.assertEqual(Skv.valueSourcesCa, {'default'; 'override'})

      % Test assertion agains reusing valueSource.
      T.verifyError(...
        @() Skv.override(111, 'default'), ...
        ?MException)
      T.verifyError(...
        @() Skv.override(222, 'override'), ...
        ?MException)
    end



    % NOTE: Operators == and ~= are not defined.
    function test_equality(T)
      Skv1 = bicas.settings.KeyValue(NaN, 'default');
      Skv2 = bicas.settings.KeyValue(NaN, 'default');
      Skv3 = bicas.settings.KeyValue(1,   'default');

      T.assertTrue( isequaln(Skv1, Skv2))
      T.assertFalse(isequaln(Skv1, Skv3))
    end



    function test0(T)

      %             % Arbitrary number output variables.
      %             function test(inputsCa, expOutputsCa)
      %                 % Pre-allocate correct size for later assignment via function.
      %                 actOutputs = cell(size(expOutputsCa));
      %
      %                 [actOutputs{:}] = FUNCTION_TO_TEST(inputsCa{:});
      %                 T.verifyEqual(actOutputs, expOutputsCa)
      %             end

      %             % One output variable.
      %             function test(inputsCa, expOutput)
      %                 actOutput = FUNCTION_TO_TEST(inputsCa{:});
      %                 T.verifyEqual(actOutput, expOutput)
      %             end

      %             function test_exc(varargin)
      %                 T.verifyError(...
      %                     @() FUNCTION_TO_TEST(varargin{:}), ...
      %                     ?MException)
      %             end
      %===================================================================


    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
