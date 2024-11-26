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

      ENV = zeros(1,0);   % Empty Logical Vector
      ELV = false(1,0);   % Empty Numeric Vector

      % Zero elements.
      testCase.test(ENV, ELV, 0, 0, ELV)
      testCase.test(ENV, ELV, 2, 0, ELV)
      testCase.test(ENV, ELV, 0, 2, ELV)

      % One element.
      testCase.test([5], [1], 0, 0, [1]);
      testCase.test([5], [1], 2, 3, [1]);
      testCase.test([5], [0], 0, 0, [0]);
      testCase.test([5], [0], 2, 3, [0]);

      % Test edges.
      testCase.test([1:13], ...
        [1,0,0,0,0,0,0,0,0,0,0,0,1], 0,   0, ...
        [1,0,0,0,0,0,0,0,0,0,0,0,1]);
      testCase.test([1:13], ...
        [1,0,0,0,0,0,0,0,0,0,0,0,1], 2,   0, ...
        [1,0,0,0,0,0,0,0,0,0,1,1,1]);
      testCase.test([1:13], ...
        [1,0,0,0,0,0,0,0,0,0,0,0,1], 0,   2, ...
        [1,1,1,0,0,0,0,0,0,0,0,0,1]);

      % "Complex" tests
      % True at edges.
      testCase.test([1:13], ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1], 0,   0, ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1]);
      testCase.test([1:13], ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1], 0.9, 0.9, ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1]);
      testCase.test([1:13], ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1], 1,   2, ...
        [1,1,1,0,1,1,1,1,1,1,0,1,1]);
      testCase.test([1:13], ...
        [1,0,0,0,0,1,0,1,0,0,0,0,1], 2,   1, ...
        [1,1,0,1,1,1,1,1,1,0,1,1,1]);

      % xMargin1/2 = Inf
      testCase.test([1:13], ...
        [0,0,0,1,0,0,0,0,0,0,0,0,0], 0,   0,   ...
        [0,0,0,1,0,0,0,0,0,0,0,0,0]);
      testCase.test([1:13], ...
        [0,0,0,1,0,0,0,0,0,0,0,0,0], 1,   inf, ...
        [0,0,1,1,1,1,1,1,1,1,1,1,1]);
      testCase.test([1:13], ...
        [0,0,0,1,0,0,0,0,0,0,0,0,0], Inf, 2, ...
        [1,1,1,1,1,1,0,0,0,0,0,0,0]);
      testCase.test([1:13], ...
        [0,0,0,1,0,0,0,0,0,0,0,0,0], Inf, inf, ...
        [1,1,1,1,1,1,1,1,1,1,1,1,1]);
    end



    function test_perm(testCase)
      % Test random x order.
      p = [3,5,9,7,4,2,1,6,13,11,12,10,8];    % Permutation
      testCase.test(...
        testCase.perm([1:13], p), ...
        testCase.perm([1,0,0,0,0,1,0,1,0,0,0,0,1], p), 1, 1, ...
        testCase.perm([1,1,0,0,1,1,1,1,1,0,0,1,1], p));
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test(testCase, x, b1, xMargin1, xMargin2, expB2)
      x     = x';
      b1    = logical(b1');
      expB2 = logical(expB2');

      actB2 = irf.utils.true_with_margin(x, b1, xMargin1, xMargin2);
      testCase.assertEqual(actB2, expB2)
    end



  end    % methods(Access=private)



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
