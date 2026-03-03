%
% matlab.unittest automatic test code for bicas.ga.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ga___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_normalize_value(T)

      function test(x, beforeValuesCa, afterValue, expX)
        actX = bicas.ga.normalize_value(x, beforeValuesCa, afterValue);
        T.assertEqual(actX, expX)
      end

      % Test the "algorithm" using numbers
      test(3, cell(0, 1), NaN, 3)
      test(3, {9}',       NaN, 3)

      test(3, {3}',     NaN, NaN)
      test(3, {2,3,4}', NaN, NaN)

      test(NaN, {NaN}', 9, 9)
      test(NaN, {3}',   9, NaN)

      % ================
      % Test non-numbers
      % ================

      % Values which commonly represent an empty GA (whether originating from
      % CDF or from BICAS code).
      GA_EMPTY_CA = {...
        [], ...
        {' '}, ...
        {'none'} ...
        }';

      test({' '},    GA_EMPTY_CA, {' '},      {' '})
      test({'02'},   GA_EMPTY_CA, {' '},      {'02'})
      test({'none'}, GA_EMPTY_CA, {' '},      {' '})
      test([],       GA_EMPTY_CA, {' '},      {' '})
      test({' '},    GA_EMPTY_CA, cell(0, 1), cell(0, 1))

      % Not normalizing (as intended!)
      test({'asd', 'none'}', GA_EMPTY_CA, {' '}, {'asd', 'none'}')
    end



    function test_normalize_empty_column_array(T)
      actGaCa = bicas.ga.normalize_empty_column_array(cell(0, 1), '-');
      T.assertEqual(actGaCa, {'-'})

      actGaCa = bicas.ga.normalize_empty_column_array({'abc'}, '-');
      T.assertEqual(actGaCa, {'abc'})

      actGaCa = bicas.ga.normalize_empty_column_array({'abc'; 'ABC'}, '-');
      T.assertEqual(actGaCa, {'abc'; 'ABC'})
    end



  end    % methods(Test)



end
