%
% matlab.unittest automatic test code for bicas.tf.apply_TF_kernel().
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-11
%
classdef apply_TF_kernel___UTEST < matlab.unittest.TestCase



  properties(TestParameter)
    % All legal values for argument edgePolicy.
    EDGE_POLICY = {'ZEROS', 'MIRROR'}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Tests which are independent of edgePolicy.
    function test_0(T, EDGE_POLICY)
      T.test( ...
        [2], [3], 1, EDGE_POLICY, ...
        [3]*2)

      T.test( ...
        [0,0,2,0,0,0], [3], 1, EDGE_POLICY, ...
        [0,0,3,0,0,0]*2)

      %===============================================
      % Even & odd length kernels.
      % Kernel origin at beginning and end of kernel.
      %===============================================
      T.test( ...
        [0,0,1,0,0,0], [1,2], 1, EDGE_POLICY, ...
        [0,0,1,2,0,0]*1)
      T.test( ...
        [0,0,1,0,0,0], [1,2], 2, EDGE_POLICY, ...
        [0,1,2,0,0,0]*1)
      T.test( ...
        [0,0,1,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
        [0,0,1,2,3,0]*1)
      T.test( ...
        [0,0,1,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
        [1,2,3,0,0,0]*1)

      %=====================================
      % Multiple non-unit, non-zero samples
      %=====================================
      T.test( ...
        [0,2,3,0,0,0], [1,2,3], 2, EDGE_POLICY, ...
        [1,2,3,0,0,0]*2 + ...
        [0,1,2,3,0,0]*3)

      % Empty kernel ==> Error
      % NOTE: Check error ID to differentiate error from any other empty
      % kernel-caused error.
      T.verifyError(...
        @() bicas.tf.apply_TF_kernel([2]', zeros(0,1), 1, EDGE_POLICY), ...
        'BICAS:Assertion:IllegalArgument')
    end



    % Special tests for edgeCase == 'ZEROS'
    %
    function test_ZEROS(T)
      y1      = [1,0,0,0,0,2];
      yKernel = [1,2,3];
      y2p     = conv(y1, yKernel);

      for iKc = 1:3
        test(T, ...
          y1, yKernel, iKc, 'ZEROS', ...
          y2p(iKc + [0:5]))
      end

      % Kernel longer than signal.
      test(T, ...
        [1,0,0], [1,2,3,4], 2, 'ZEROS', ...
        [2,3,4])
    end



    % Special tests for edgeCase == 'CYCLIC'
    %
    function test_CYCLIC(T)
      for iKo = 1:3
        T.test_pad( ...
          [0,0,2], [1,0,0,0,0,2], [1,0,0], [1,2,3], iKo, 'CYCLIC')
      end

      T.test_pad( ...
        [3,4], [2,1,0,0,3,4], [2,1], [1,2,3,4,5], 3, 'CYCLIC')

      % Kernel longer than signal (but still legal).
      T.test_pad( ...
        [1,0,0], [1,0,0], [1,0,0], [1,2,3,4], 3, 'CYCLIC')
    end



    % Special tests for edgeCase == 'MIRROR'
    %
    function test_MIRROR(T)
      for iKo = 1:3
        T.test_pad( ...
          [0,0,1], [1,0,0,0,0,2], [2,0,0], [1,2,3], iKo, 'MIRROR')
      end

      T.test_pad( ...
        [1,2], [2,1,0,0,3,4], [4,3], [1,2,3,4,5], 3, 'MIRROR')

      % Kernel longer than signal (but still legal).
      T.test_pad( ...
        [0,0,1], [1,0,0], [0,0,1], [1,2,3,4], 2, 'MIRROR')
    end



    function test_NaN(T, EDGE_POLICY)
      N = NaN;   % Shorthand.

      %===============
      % NaN in signal
      %===============
      T.test( ...
        [0,0,0,N,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
        [0,N,N,N,0,0,0])
      T.test( ...
        [0,0,0,N,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
        [0,0,0,N,N,N,0])

      %===============
      % NaN in kernel
      %===============
      T.test( ...
        [0,0,0,0,0,0,0], [N,0,0], 3, EDGE_POLICY, ...
        [N,N,N,N,N,N,N])
      T.test( ...
        [0,0,0,0,0,0,0], [N,0,0], 1, EDGE_POLICY, ...
        [N,N,N,N,N,N,N])
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test(T, y1, yKernel, iKernelOrigin, edgePolicy, expY2)

      % NOTE: Transposes vectors.
        actY2 = bicas.tf.apply_TF_kernel(...
        y1', yKernel', iKernelOrigin, edgePolicy);

      T.verifyEqual(actY2, expY2')
    end



    % Manually pad y1.
    % NOTE: Function derives the expected output itself, using conv().
    %
    function test_pad(T, y1a, y1, y1b, yKernel, iKernelOrigin, edgePolicy)

      assert(isrow(y1a))
      assert(isrow(y1))
      assert(isrow(y1b))
      assert(isrow(yKernel))

      y1a     = y1a';
      y1      = y1';
      y1b     = y1b';
      yKernel = yKernel';

      y1p   = [y1a; y1; y1b];
      % NOTE: Using conv() when applying the kernel.
      y2p   = conv(y1p, yKernel);

      iBegin = length(y1a) + iKernelOrigin;
      iEnd   = iBegin + length(y1) - 1;
      y2exp = y2p(iBegin : iEnd);

      T.test(...
        y1', yKernel', iKernelOrigin, edgePolicy, ...
        y2exp')
    end



  end    % methods(Static, Access=private)



end
