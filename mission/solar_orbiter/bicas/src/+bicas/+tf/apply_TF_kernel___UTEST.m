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
        EDGE_POLICY = {'zeros', 'mirror'}
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase, EDGE_POLICY)
            import bicas.tf.apply_TF_kernel___UTEST.test

            test(testCase, ...
                [2], [3], 1, EDGE_POLICY, ...
                [3]*2)

            test(testCase, ...
                [0,0,2,0,0,0], [3], 1, EDGE_POLICY, ...
                [0,0,3,0,0,0]*2)

            %===============================================
            % Even & odd length kernels.
            % Kernel origin at beginning and end of kernel.
            %===============================================
            test(testCase, ...
                [0,0,1,0,0,0], [1,2], 1, EDGE_POLICY, ...
                [0,0,1,2,0,0]*1)
            test(testCase, ...
                [0,0,1,0,0,0], [1,2], 2, EDGE_POLICY, ...
                [0,1,2,0,0,0]*1)
            test(testCase, ...
                [0,0,1,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
                [0,0,1,2,3,0]*1)
            test(testCase, ...
                [0,0,1,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
                [1,2,3,0,0,0]*1)

            %=====================================
            % Multiple non-unit, non-zero samples
            %=====================================
            test(testCase, ...
                [0,2,3,0,0,0], [1,2,3], 2, EDGE_POLICY, ...
                [1,2,3,0,0,0]*2 + ...
                [0,1,2,3,0,0]*3)

            % Empty kernel ==> Error
            % NOTE: Check error ID to differentiate error from any other empty
            % kernel-caused error.
            testCase.verifyError(...
                @() bicas.tf.apply_TF_kernel([2]', zeros(0,1), 1, EDGE_POLICY), ...
                'BICAS:Assertion:IllegalArgument')
        end



        % Special tests for edgeCase == 'zeros'
        %
        function test_zeros(testCase)
            import bicas.tf.apply_TF_kernel___UTEST.test

            y1      = [1,0,0,0,0,2];
            yKernel = [1,2,3];
            y2p     = conv(y1, yKernel);

            for iKc = 1:3
                test(testCase, ...
                    y1, yKernel, iKc, 'zeros', ...
                    y2p(iKc + [0:5]))
            end

            % Kernel longer than signal.
            test(testCase, ...
                [1,0,0], [1,2,3,4], 2, 'zeros', ...
                [2,3,4])
        end



        % Special tests for edgeCase == 'cyclic'
        %
        function test_cyclic(testCase)
            import bicas.tf.apply_TF_kernel___UTEST.test
            import bicas.tf.apply_TF_kernel___UTEST.test_pad

            for iKo = 1:3
                test_pad(testCase, ...
                    [0,0,2], [1,0,0,0,0,2], [1,0,0], [1,2,3], iKo, 'cyclic')
            end

            test_pad(testCase, ...
                [3,4], [2,1,0,0,3,4], [2,1], [1,2,3,4,5], 3, 'cyclic')

            % Kernel longer than signal (but still legal).
            test_pad(testCase, ...
                [1,0,0], [1,0,0], [1,0,0], [1,2,3,4], 3, 'cyclic')
        end



        % Special tests for edgeCase == 'mirror'
        %
        function test_mirror(testCase)
            import bicas.tf.apply_TF_kernel___UTEST.test
            import bicas.tf.apply_TF_kernel___UTEST.test_pad

            for iKo = 1:3
                test_pad(testCase, ...
                    [0,0,1], [1,0,0,0,0,2], [2,0,0], [1,2,3], iKo, 'mirror')
            end

            test_pad(testCase, ...
                [1,2], [2,1,0,0,3,4], [4,3], [1,2,3,4,5], 3, 'mirror')

            % Kernel longer than signal (but still legal).
            test_pad(testCase, ...
                [0,0,1], [1,0,0], [0,0,1], [1,2,3,4], 2, 'mirror')
        end



        function test_NaN(testCase, EDGE_POLICY)
            import bicas.tf.apply_TF_kernel___UTEST.test
            N = NaN;   % Shorthand.

            %===============
            % NaN in signal
            %===============
            test(testCase, ...
                [0,0,0,N,0,0,0], [1,2,3], 3, EDGE_POLICY, ...
                [0,N,N,N,0,0,0])
            test(testCase, ...
                [0,0,0,N,0,0,0], [1,2,3], 1, EDGE_POLICY, ...
                [0,0,0,N,N,N,0])

            %===============
            % NaN in kernel
            %===============
            test(testCase, ...
                [0,0,0,0,0,0,0], [N,0,0], 3, EDGE_POLICY, ...
                [N,N,N,N,N,N,N])
            test(testCase, ...
                [0,0,0,0,0,0,0], [N,0,0], 1, EDGE_POLICY, ...
                [N,N,N,N,N,N,N])
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
    %methods(Static, Access=public)



        function test(testCase, ...
                y1, yKernel, iKernelOrigin, edgePolicy, ...
                expOutput)

            % NOTE: Transposes vectors.
            testCase.verifyEqual(...
                bicas.tf.apply_TF_kernel(...
                y1', yKernel', iKernelOrigin, edgePolicy),...
                expOutput')
        end



        % Manually pad y1.
        % NOTE: Function derives the expected output itself, using conv().
        %
        function test_pad(testCase, ...
                y1a, y1, y1b, yKernel, iKernelOrigin, edgePolicy)

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

            testCase.verifyEqual(...
                bicas.tf.apply_TF_kernel(...
                    y1, yKernel, iKernelOrigin, edgePolicy),...
                y2exp)
        end



    end    % methods(Static, Access=private)



end
