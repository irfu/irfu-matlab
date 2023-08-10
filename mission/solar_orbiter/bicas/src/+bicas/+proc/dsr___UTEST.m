%
% matlab.unittest automatic test code for bicas.proc.dsr.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-10 from older test code.
%
classdef dsr___UTEST < matlab.unittest.TestCase



    properties(Constant)
        L = bicas.Logger('none', false);
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test_downsample_sci_ZV(testCase)
            
            % function [zvMed, zvMstd] = downsample_sci_ZV(...
            %         zv, nMinReqRecords, iRecordsInBinCa, L)

            % Arbitrary number output variables.
            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));
                
                [actOutputs{:}] = bicas.proc.dsr.downsample_sci_ZV(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            
            % Create test with exactly ONE BIN.
            function test_1_bin(zv, nMinReqSamples, med, mstd)
                test(...
                    {zv, nMinReqSamples, {1:size(zv,1)}, testCase.L}, ...
                    {med, mstd});
            end

            % Create test with N BINS (i.e. an arbitrary call).
            function test_N_bins(zv, nMinReqSamples, iRecordsDsrCa, med, mstd)
                assert(isrow(iRecordsDsrCa))
                bicas.proc.dsr___UTEST.assert_iRecordsDsrCa(iRecordsDsrCa, zv)

                test({zv, nMinReqSamples, iRecordsDsrCa', testCase.L}, ...
                    {med, mstd});
            end
            %===================================================================

            AS10 = zeros(1,0);   % AS10 = Array Size 1x0
            N    = NaN;


            % Empty data, zero records, zero columns.
            test_1_bin(zeros(0,0), 0, AS10, AS10);
            test_1_bin(zeros(2,0), 0, AS10, AS10);
            test_1_bin(zeros(0,2), 0, NaN(1,2), NaN(1,2));

            % mstd=0
            test_1_bin([1,2,3              ], 0, [1,2,3], [nan,nan,nan]);
            test_1_bin([1,2,3; 1,2,3       ], 0, [1,2,3], [0,0,0]);
            test_1_bin([1,2,3; 1,2,3; 1,2,3], 0, [1,2,3], [0,0,0]);
            test_1_bin([1    ; 1    ; 1    ], 0, [1],     [0]);

            % Test nMinReqSamples
            % -------------------
            % Enough samples
            test_1_bin([1,2,3; 1,2,3; 1,2,3], 3, [1,2,3], [0,0,0]);
            % Not enough samples
            test_1_bin([1,2,3; 1,2,3; 1,2,3], 4, [N,N,N], [N,N,N]);
            % Not enough samples if removes NaN.
            test_1_bin([1,2,3; 1,2,3; 1,2,N],        3, [1,2,N], [0,0,N]);
            % Enough samples even if removes NaN.
            test_1_bin([1,2,3; 1,2,3; 1,2,N; 1,2,3], 3, [1,2,3], [0,0,0]);



            % Average of two values (special case)
            test_1_bin([1,2,3; 2,3,4], 0, [1.5, 2.5, 3.5], sqrt(0.5)*[1,1,1]);
            % Nominal median
            test_1_bin([1;2;10],       0, [2], sqrt( (1^2+0^2+8^2)/2 ));



            test_N_bins([1,2; 2,3], 1, {1, 2}, [1,2; 2,3], NaN(2,2))

            test_N_bins(...
                [1, 2; ...
                 2, 3; ...
                 3, 4; ...
                 4, 5; ...
                 5, 6], 1, {1:2, 3:5}, ...
                [1.5, 2.5; ...
                 4,   5   ], ...
                [sqrt(2*0.5^2)/1 * [1,1]; ...
                 sqrt(4*1^2  )/2 * [1,1]])

            test_N_bins(...
                [1, 2; ...
                 2, 3; ...
                 3, 4; ...
                 4, 5; ...
                 5, 6], 1, {1:2, [], 3:5}, ...
                [1.5, 2.5; ...
                 NaN, NaN; ...
                 4,   5], ...
                [sqrt(2*0.5^2)/1 * [1,1]; ...
                 NaN, NaN; 
                 sqrt(4*1^2  )/2 * [1,1]])

            % Higher threshold
            test_N_bins([1,2; 2,3; 3,4; 4,5; 5,6], 3, {1:2, [], 3:5}, ...
                [NaN, NaN; ...
                 NaN, NaN; ...
                 4, 5], ...
                [NaN, NaN; ...
                 NaN, NaN; 
                 sqrt(4*1^2  )/2 * [1,1]])

        end
        
        
        
        function test_get_downsampling_bins(testCase)
            
            function add_test(...
                    zvAllUtcCa, ...
                    boundaryRefUtc, ...
                    intervalLengthWolsNs, ...
                    timestampPosWolsNs, ...
                    ...
                    zvIntervalsUtcCa, ...
                    iRecordsCa, ...
                    binSizeArrayNs)

                inputsCa = {...
                    spdfparsett2000(zvAllUtcCa), ...
                    spdfparsett2000(boundaryRefUtc), ...
                    int64(intervalLengthWolsNs), ...
                    int64(timestampPosWolsNs), ...
                    testCase.L};
                expOutputsCa = {...
                    spdfparsett2000(zvIntervalsUtcCa), ...
                    iRecordsCa(:), ...
                    int64(binSizeArrayNs(:))};
                
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));
                
                [actOutputs{:}] = bicas.proc.dsr.get_downsampling_bins(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
                
            end
            %===================================================================

            % ECA = Empty Column (size=0) Array
            ECA = zeros(0,1);

            % TEST: No timestamps/empty Epoch.
            add_test({}, '2020-01-01T00:00:00', 10e9, 5e9, {}, {}, ECA);

            % TEST: One timestamps.
            add_test({'2020-01-01T00:00:06'}, '2020-01-01T00:00:00', 10e9, 5e9, {'2020-01-01T00:00:05'}, {1}, [10e9]);
            add_test({'2020-01-01T00:00:06'}, '2020-01-01T00:00:01', 10e9, 3e9, {'2020-01-01T00:00:04'}, {1}, [10e9]);

            % TEST: Varying number of timestamps in each interval
            add_test(...
                {...
                '2020-01-01T00:01:06', ...
                '2020-01-01T00:01:12', ...
                '2020-01-01T00:01:18'
                }, ...
                '2020-01-01T00:00:00', 10e9, 5e9, ...
                {...
                '2020-01-01T00:01:05', ...
                '2020-01-01T00:01:15'...
                }, {1, [2,3]'}, [10e9, 10e9]);

            % TEST: Total interval over LEAP SECOND (end of 2016).
            % TEST: Data gap (two intervals without timestamps).
            % TEST: Boundary reference timestamp far from data.
            % TEST: Generally convoluted case...
            add_test(...
                {...
                '2016-12-31T23:59:46', ...
                '2016-12-31T23:59:51', ...
                '2016-12-31T23:59:59', ...
                '2017-01-01T00:00:02', ...
                '2017-01-01T00:00:03', ...
                '2017-01-01T00:00:04', ...
                '2017-01-01T00:00:34', ...
                }, ...
                '2020-01-01T00:00:05', 10e9, 5e9, ...
                {...
                '2016-12-31T23:59:50', ...
                '2017-01-01T00:00:00'...
                '2017-01-01T00:00:10', ...
                '2017-01-01T00:00:20'...
                '2017-01-01T00:00:30'...
                }, {[1,2]', [3,4,5,6]', ECA, ECA, [7]}, [10e9, 11e9, 10e9, 10e9, 10e9]);

        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Internal utility function
        function assert_iRecordsDsrCa(iRecordsDsrCa, zv)

            % ~Normalize
            iRecordsDsrNormCa = arrayfun(...
                @(s) (s{1}(:)'), ...
                iRecordsDsrCa, ...
                'UniformOutput', false);

            iRecordsArray = [iRecordsDsrNormCa{:}];

            assert(all(unique(iRecordsArray) == 1:size(zv, 1)))
        end
        
        
        
    end    % methods(Static, Access=private)

    
    
end
