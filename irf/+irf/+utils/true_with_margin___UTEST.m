%
% matlab.unittest automatic test code for irf.utils.true_with_margin().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef true_with_margin___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Iterate over tests using TestParameter.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      % PROPOSAL: Test more unsorted x.

      function test(x, b1, xMargin, expB2)
        x     = x';
        b1    = logical(b1');
        expB2 = logical(expB2');

        actB2 = irf.utils.true_with_margin(x, b1, xMargin);
        testCase.assertEqual(actB2, expB2)
      end

      ENV = zeros(1,0);   % Empty Logical Vector
      ELV = false(1,0);   % Empty Numeric Vector

      test(ENV, ELV, 0, ELV)
      test([5], [1], 0, [1]);
      test([5], [1], 2, [1]);
      test([5], [0], 0, [0]);
      test([5], [0], 2, [0]);

      % Test different xMargin.
      % Test true at edges.
      test([1:13], [1,0,0,0,0,1,0,1,0,0,0,0,1], 0,   [1,0,0,0,0,1,0,1,0,0,0,0,1]);
      test([1:13], [1,0,0,0,0,1,0,1,0,0,0,0,1], 0.9, [1,0,0,0,0,1,0,1,0,0,0,0,1]);
      test([1:13], [1,0,0,0,0,1,0,1,0,0,0,0,1], 1,   [1,1,0,0,1,1,1,1,1,0,0,1,1]);
      test([1:13], [1,0,0,0,0,1,0,1,0,0,0,0,1], 2,   [1,1,1,1,1,1,1,1,1,1,1,1,1]);

      test([1:13], [1,0,0,0,0,0,0,0,0,0,0,0,1], 0,   [1,0,0,0,0,0,0,0,0,0,0,0,1]);
      test([1:13], [1,0,0,0,0,0,0,0,0,0,0,0,1], 1,   [1,1,0,0,0,0,0,0,0,0,0,1,1]);
      test([1:13], [1,0,0,0,0,0,0,0,0,0,0,0,1], 2,   [1,1,1,0,0,0,0,0,0,0,1,1,1]);

      test([1:13], [0,0,0,1,0,0,0,0,0,0,0,0,1], Inf, [1,1,1,1,1,1,1,1,1,1,1,1,1]);

      % Test random x order.
      p = [3,5,9,7,4,2,1,6,13,11,12,10,8];    % Permutation
      test(...
        testCase.perm([1:13], p), ...
        testCase.perm([1,0,0,0,0,1,0,1,0,0,0,0,1],p), 1, ...
        testCase.perm([1,1,0,0,1,1,1,1,1,0,0,1,1],p));
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function x = perm(x, p)
      % Convert to static function.

      % Assert that p=permutation.
      assert(isequal(sort(p), 1:numel(p)))

      x = x(p);
    end



  end    % methods(Static, Access=private)



end
