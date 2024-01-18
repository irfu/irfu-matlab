%
% matlab.unittest automatic test code for bicas.proc.L1L2.cal.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16
%
classdef utils___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_get_calibration_time(testCase)

            function test(Epoch, CalibEpochList, expOutput)
                % NOTE: Converting arguments to int64, transposing.
                actOutput = bicas.proc.L1L2.cal.utils.get_calibration_time(...
                        int64(Epoch)', ...
                        int64(CalibEpochList)');

                testCase.verifyEqual(...
                    actOutput, ...
                    expOutput')
            end

            function test_exc(Epoch, CalibEpochList)
                f = @() bicas.proc.L1L2.cal.utils.get_calibration_time(...
                        int64(Epoch)', int64(CalibEpochList)');

                testCase.verifyError(...
                    f, ?MException);
            end

            %===================================================================

            EV = zeros(1,0);

            test(EV,     [3],   EV)
            test(EV,     [3,7], EV)

            test([5:10], [3]   , [1,1,1,1,1,1]);
            test([5:10], [3, 7], [1,1,2,2,2,2]);

            % CalibEpochList must not be empty.
            test_exc([1:10], []);
            % Requesting indices before first timestamp.
            test_exc([1:10], [3])
            test_exc([1:10], [3, 7]);

        end



    end    % methods(Test)



end
