%
% UNFINISHED / EXPERIMENTAL
%
% matlab.unittest automatic test code for irf.utils.mixed_radix().
%
% NOTE: Hard-coded "mrd" and "bases" arguments  are on more user-friendly
% format, and not the format used for arguments to the functions being tested.
% iDigit=1 <=> Most significant digit/base.
% base : row vector
%
% NOTE: Test functions always have arguments in order (n, mrd, bases) regardless
% of direction of conversion.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-18, using code from other test code.
%
classdef mixed_radix___UTEST < matlab.unittest.TestCase



  properties(TestParameter)
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_I2M(testCase)

      function test(n, expMrd, bases)
        expMrd = int64(fliplr(expMrd));
        bases  = fliplr(bases)';

        actMrd = irf.utils.mixed_radix.integer_to_MRD(...
          n, bases);
        testCase.verifyEqual(actMrd, expMrd)
      end

      function test_exc(varargin)
        testCase.verifyError(...
          @() irf.utils.mixed_radix.integer_to_MRD(varargin{:}), ...
          ?MException)
      end
      %===================================================================

      % Negative integers.
      % NOTE: Conversion (integer <--> MRD) can not be symmetric for
      % negative integers since MRD always converts to a non-negative
      % integer.
      test( -1, [0,1,2], [1,2,3]);
      test(  6, [0,0,0], [1,2,3]);
      test( -6, [0,0,0], [1,2,3]);
      test( 10, [0,1,1], [1,2,3]);
      test(-10, [0,0,2], [1,2,3]);



      % Check that submitting double-typed arguments is OK.
      actMrd = irf.utils.mixed_radix.integer_to_MRD(...
        3, fliplr([2,2,2])');
      testCase.verifyEqual(actMrd, int64(fliplr([0,1,1])))

      % Base<=0 illegal.
      test_exc(5, [2, 0,2]')
      test_exc(5, [2,-2,2]')
    end



    function test_M2I(testCase)

      function test(n, mrd, bases)
        n     = int64(n);
        mrd   = int64(fliplr(mrd));
        bases = fliplr(bases)';

        actN = irf.utils.mixed_radix.MRD_to_integer(...
          mrd, bases);
        testCase.verifyEqual(actN, n)

      end

      function test_exc(varargin)
        testCase.verifyError(...
          @() irf.utils.mixed_radix.MRD_to_integer(varargin{:}), ...
          ?MException)
      end
      %===================================================================

      % Check that submitting double-typed arguments is OK.
      actN = irf.utils.mixed_radix.MRD_to_integer(...
        fliplr([1,1,0]), fliplr([2,2,2])');
      testCase.verifyEqual(actN, int64(6))

      % Base<=0 illegal.
      test_exc([0, 0, 0], [2, 0,2]')
      test_exc([0, 0, 0], [2,-2,2]')

      % MRD >= base.
      test(10, [0, 4, 2], [4,3,2])
    end



    % Check that conversions in both directions are symmetric.
    function test_I2M_M2I(testCase)

      function test(n, mrd, bases)
        n     = int64(n);
        mrd   = int64(fliplr(mrd));
        bases = fliplr(bases)';

        actMrd = irf.utils.mixed_radix.integer_to_MRD(...
          n, bases);
        testCase.verifyEqual(actMrd, mrd)


        actN = irf.utils.mixed_radix.MRD_to_integer(...
          mrd, bases);
        testCase.verifyEqual(actN, n)

      end

      %===================================================================

      % Array size edge cases.
      test(zeros(0,1), zeros(0,0), zeros(1,0));   % Zero bases, zero numbers
      test(zeros(3,1), zeros(3,0), zeros(1,0));   % Zero bases
      test(zeros(0,1), zeros(0,3), [1,2,3]);      % Zero numbers

      % Base edge cases.
      test([0],   [0,0,0,0], [1,1,1,1]);
      test([2+1], [0,1,1,0], [1,2,2,1]);

      % Regular conversions
      test(5, [0,1,0,1], [2,2,2,2]);
      test([[1;3]*30+[2;2]*6+[3;1]], [1,2,3; 3,2,1], [4,5,6]);
      test([8+2+1],       [1,0,1,1], [2,2,2,2]);
      test([4],           [0,1,0,0], [2,2,2,2]);
      test([3*6+2*2+1*1], [3,2,1],   [4,3,2]);
      test([3;1]*5*6 + [2;2]*6 + [1;3], [3,2,1; 1,2,3], [4,5,6]);

      % Test multiple mRDs/integers.
      test([3;8], [0,1,1; 1,1,0], [4,3,2])
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
