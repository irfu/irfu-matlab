%
% matlab.unittest automatic test code for bicas.tf.kernel.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-11
%
classdef kernel___UTEST < matlab.unittest.TestCase



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



    % Tests which are independent of edgePolicy by their nature.
    function test_apply_TF___0(T, EDGE_POLICY)
      T.test_apply_TF( ...
        [2], [3], 1, EDGE_POLICY, ...
        [3]*2)

      T.test_apply_TF( ...
        [0,0,2,0,0,0], [3], 1, EDGE_POLICY, ...
        [0,0,3,0,0,0]*2)

      %===============================================
      % Even & odd length kernels.
      % Kernel origin at beginning and end of kernel.
      %===============================================
      T.test_apply_TF( ...
        [0,0,1,0,0,0], [1,2], 1, EDGE_POLICY, ...
        [0,0,1,2,0,0]*1)
      T.test_apply_TF( ...
        [0,0,1,0,0,0], [1,2], 2, EDGE_POLICY, ...
        [0,1,2,0,0,0]*1)
      T.test_apply_TF( ...
        [0,0,1,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
        [0,0,1,2,3,0]*1)
      T.test_apply_TF( ...
        [0,0,1,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
        [1,2,3,0,0,0]*1)

      %=====================================
      % Multiple non-unit, non-zero samples
      %=====================================
      T.test_apply_TF( ...
        [0,2,3,0,0,0], [1,2,3], 2, EDGE_POLICY, ...
        [1,2,3,0,0,0]*2 + ...
        [0,1,2,3,0,0]*3)

      % Empty kernel ==> Error
      % NOTE: Check error ID to differentiate error from any other empty
      % kernel-caused error.
      T.assertError(...
        @() bicas.tf.kernel.apply_kernel([2]', zeros(0,1), 1, EDGE_POLICY), ...
        'BICAS:Assertion:IllegalArgument')
    end



    % Special tests for edgeCase == 'ZEROS'
    %
    function test_apply_TF___ZEROS(T)
      y1      = [1,0,0,0,0,2];
      yKernel = [1,2,3];
      y2p     = conv(y1, yKernel);

      for iKc = 1:3
        T.test_apply_TF( ...
          y1, yKernel, iKc, 'ZEROS', ...
          y2p(iKc + [0:5]))
      end

      % Kernel has same length as signal (max length).
      T.test_apply_TF( ...
        [1,0,0], [1,2,3], 2, 'ZEROS', ...
        [2,3,0])
    end



    % Special tests for edgeCase == 'CYCLIC'
    %
    function test_apply_TF___CYCLIC(T)
      for iKo = 1:3
        T.test_apply_TF_pad( ...
          [0,0,2], [1,0,0,0,0,2], [1,0,0], [1,2,3], iKo, 'CYCLIC')
      end

      T.test_apply_TF_pad( ...
        [3,4], [2,1,0,0,3,4], [2,1], [1,2,3,4,5], 3, 'CYCLIC')

      % Kernel has same length as signal (max length).
      T.test_apply_TF_pad( ...
        [1,0,0], [1,0,0], [1,0,0], [1,2,3], 3, 'CYCLIC')
    end



    % Special tests for edgeCase == 'MIRROR'
    %
    function test_apply_TF___MIRROR(T)
      for iKo = 1:3
        T.test_apply_TF_pad( ...
          [0,0,1], [1,0,0,0,0,2], [2,0,0], [1,2,3], iKo, 'MIRROR')
      end

      T.test_apply_TF_pad( ...
        [1,2], [2,1,0,0,3,4], [4,3], [1,2,3,4,5], 3, 'MIRROR')

      % Kernel has same length as signal (max length).
      T.test_apply_TF_pad( ...
        [0,0,1], [1,0,0], [0,0,1], [1,2,3], 2, 'MIRROR')
    end



    function test_apply_TF___NaN(T, EDGE_POLICY)
      N = NaN;   % Shorthand.

      %===============
      % NaN in signal
      %===============
      T.test_apply_TF( ...
        [0,0,0,N,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
        [0,N,N,N,0,0,0])
      T.test_apply_TF( ...
        [0,0,0,N,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
        [0,0,0,N,N,N,0])

      %===============
      % NaN in kernel
      %===============
      T.test_apply_TF( ...
        [0,0,0,0,0,0,0], [N,0,0], 3, EDGE_POLICY, ...
        [N,N,N,N,N,N,N])
      T.test_apply_TF( ...
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



    function test_apply_TF(T, y1, yKernel, iKernelOrigin, edgePolicy, expY2)

      % NOTE: Transposes vectors.
      actY2 = bicas.tf.kernel.apply_kernel(...
        y1', yKernel', iKernelOrigin, edgePolicy);

      T.assertEqual(actY2, expY2')
    end



    % Helper function for testing edgePolicy<>ZEROS.
    % (1) Manually pad y1 via arguments, corresponding to some edgePolicy value.
    % (2) derive expected y2 using bicas.tf.kernel.apply_kernel() with
    %     edgePolicy=ZEROS.
    % (3) Check consistency of result (2) with bicas.tf.kernel.apply_kernel()
    %     WITHOUT manual padding, but specified edgePolicy<>ZEROS.
    %
    function test_apply_TF_pad(T, ...
        y1PadA, y1, y1PadB, yKernel, iKernelOrigin, edgePolicy)
      % PROBLEM: Overlaps too much with the tested function itself?!
      %   PROPOSAL: Tests should hardcode explicit outputs?!
      %     CON: Not only, due to there being too many cases (iKernelorigin).
      %          Consistency checks between different types of calls are still
      %          useful.
      %       CON-PROPOSAL: There should be some non-ZEROS test with explicitly
      %                     hardcoded output.

      assert(isrow(y1PadA))
      assert(isrow(y1))
      assert(isrow(y1PadB))
      assert(isrow(yKernel))

      y1PadA  = y1PadA';
      y1      = y1';
      y1PadB  = y1PadB';
      yKernel = yKernel';

      y1padded = [y1PadA; y1; y1PadB];

      iBegin = length(y1PadA) + 1;
      iEnd   = iBegin + length(y1) - 1;

      % Derive expected y2 using bicas.tf.kernel.apply_kernel() using MANUALLY
      % PADDED y1.
      expY2padded = bicas.tf.kernel.apply_kernel(...
        y1padded, yKernel, iKernelOrigin, 'ZEROS');
      expY2 = expY2padded(iBegin:iEnd);

      % Call bicas.tf.kernel.apply_kernel() WITHOUT manual padding.
      T.test_apply_TF(...
        y1', yKernel', iKernelOrigin, edgePolicy, ...
        expY2')
    end



  end    % methods(Static, Access=private)



end
