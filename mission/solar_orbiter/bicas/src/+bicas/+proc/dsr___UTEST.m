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



        function test_downsample_sci_ZV(testCase)

            function test(...
                    OsrAr, nMinNfpSamplesPerBin, iRecordsInBinCa, L, ...
                    expMedianDsrAr, expMstdDsrAr)

                % NOTE: bicas.proc.dsr.downsample_sci_ZV() only accepts
                % double-typed FPAs. Therefore only needs to test such.
                OsrFpa          = bicas.utils.FPArray(OsrAr,          'FILL_VALUE', NaN);
                ExpMedianDsrFpa = bicas.utils.FPArray(expMedianDsrAr, 'FILL_VALUE', NaN);
                ExpMstdDsrFpa   = bicas.utils.FPArray(expMstdDsrAr,   'FILL_VALUE', NaN);

                [ActMedianDsrFpa, ActMstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                    OsrFpa, nMinNfpSamplesPerBin, iRecordsInBinCa, L);

                testCase.assertEqual(ActMedianDsrFpa, ExpMedianDsrFpa)
                testCase.assertEqual(ActMstdDsrFpa,   ExpMstdDsrFpa)
            end



            % Create test with exactly ONE BIN.
            function test_1_bin(osrAr, nMinReqSamples, expMedAr, expMstdAr)
                test(...
                    osrAr, nMinReqSamples, {1:size(osrAr,1)}, testCase.L, ...
                    expMedAr, expMstdAr);
            end



            % Create test with N BINS (i.e. an arbitrary call).
            function test_N_bins(zv, nMinReqSamples, iRecordsDsrCa, med, mstd)
                assert(isrow(iRecordsDsrCa))
                bicas.proc.dsr___UTEST.assert_iRecordsDsrCa(iRecordsDsrCa, zv)

                test(...
                    zv, nMinReqSamples, iRecordsDsrCa', testCase.L, ...
                    med, mstd);
            end



            AS10 = zeros(1,0);   % AS10 = Array Size 1x0
            N    = NaN;

            %====================
            % Tests with one bin
            %====================
            % Empty data, zero records, zero columns.
            test_1_bin(zeros(0,0), 0, AS10, AS10);
            test_1_bin(zeros(2,0), 0, AS10, AS10);
            test_1_bin(zeros(0,2), 0, NaN(1,2), NaN(1,2));

            % mstd=0
            test_1_bin([1,2,3              ], 0, [1,2,3], [NaN,NaN,NaN]);
            test_1_bin([1,2,3; 1,2,3       ], 0, [1,2,3], [0,0,0]);
            test_1_bin([1,2,3; 1,2,3; 1,2,3], 0, [1,2,3], [0,0,0]);
            test_1_bin([1    ; 1    ; 1    ], 0, [1],     [0]);

            % Test nMinReqSamples
            % -------------------
            % Enough samples
            test_1_bin([1,2,3; 1,2,3; 1,2,3],        3, [1,2,3], [0,0,0]);
            % Not enough samples
            test_1_bin([1,2,3; 1,2,3; 1,2,3],        4, [N,N,N], [N,N,N]);
            % Not enough samples if removes NaN.
            test_1_bin([1,2,3; 1,2,3; 1,2,N],        3, [1,2,N], [0,0,N]);
            % Enough samples even if removes NaN.
            test_1_bin([1,2,3; 1,2,3; 1,2,N; 1,2,3], 3, [1,2,3], [0,0,0]);

            % Median of two values ==> Average (special case)
            test_1_bin([1,2,3; 2,3,4], 0, [1.5, 2.5, 3.5], sqrt(0.5)*[1,1,1]);
            % Nominal median
            test_1_bin([1;2;10],       0, [2], sqrt( (1^2+0^2+8^2)/2 ));



            %===================
            % Tests with N bins
            %===================
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
            test_N_bins(...
                [1,2;
                 2,3;
                 3,4;
                 4,5;
                 5,6], 3, {1:2, [], 3:5}, ...
                [NaN, NaN; ...
                 NaN, NaN; ...
                 4, 5], ...
                [NaN, NaN; ...
                 NaN, NaN;
                 [1,1] * sqrt(4*1^2  )/2])
        end



        function test_downsample_ZV_minimum(testCase)

            % Function handle to function to be tested, so that one can easiy
            % switch to other implementations of the same function for testing,
            % e.g. for testing performance for different implementations
            FH_CA = {};
%             FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_minimum_W_FPAs;
%             FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_minimum_W_INNER_ARRAYS;
            FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_minimum;

            function test(inAr, inFv, iRecordsInBinCa, expDsrAr, expDsrFv)
                Fpa       = bicas.utils.FPArray(inAr,     'FILL_VALUE', inFv);
                ExpDsrFpa = bicas.utils.FPArray(expDsrAr, 'FILL_VALUE', expDsrFv);
                ActDsrFpa = fh(Fpa, iRecordsInBinCa);

                testCase.assertEqual(ActDsrFpa, ExpDsrFpa)
            end

            for iFh = 1:numel(FH_CA)
                fh = FH_CA{iFh};

                % Empty
                test(ones(0, 1), NaN,        cell(0, 1), ones(0, 1), NaN)
                % One sample, non-double type
                test(uint8(3),   uint8(255), {1},        uint8(3),   uint8(255))
                % Three, different sized bins
                test([3;4;5; 2;3; 9],   NaN,  {1:3; 4:5; 6},    [3;2;9],  NaN)

                test([NaN],     NaN,  {1},      [NaN],  NaN)
                test([2;NaN],   NaN,  {1:2},    [2],    NaN)

                test([2;NaN; NaN],   NaN,  {1:2; 3},    [2; NaN],  NaN)

                % Size-zero bin.
                test(ones(0, 1),   NaN,  {1:0},         [NaN],         NaN)
                test([2;3; NaN],   NaN,  {1:2; 1:0; 3}, [2; NaN; NaN], NaN)

            end
        end



        function downsample_ZV_bitmask(testCase)
            % Function handle to function to be tested, so that one can easiy
            % switch to other implementations of the same function for testing.
            FH_CA = {};
            %FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_bitmask_W_FPAs;
            %FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_bitmask_W_INNER_ARRAYS;
            FH_CA{end+1} = @bicas.proc.dsr.downsample_ZV_bitmask;

            function test(inAr, iRecordsInBinCa, expDsrAr)
                inAr     = uint8(inAr);
                inFv     = uint8(255);
                expDsrAr = uint8(expDsrAr);
                expDsrFv = uint8(255);

                Fpa       = bicas.utils.FPArray(inAr,     'FILL_VALUE', inFv);
                ExpDsrFpa = bicas.utils.FPArray(expDsrAr, 'FILL_VALUE', expDsrFv);

                ActDsrFpa = fh(Fpa, iRecordsInBinCa);

                testCase.assertEqual(ActDsrFpa, ExpDsrFpa)
            end

            for iFh = 1:numel(FH_CA)
                fh = FH_CA{iFh};

                % Zero samples, zero bins
                test(ones(0, 1), cell(0, 1), ones(0, 1))
                % One sample
                test(3,          {1},        3)

                % One bin, multiple samples
                test([1;2;4],    {1:3},      7)

                % Three, different sized bins
                test([4;5;6; 2;3; 9],  {1:3; 4:5; 6},    [7;3;9])

                test([255],     {1},      [255])
                test([2;255],   {1:2},    [2])

                test([2;255; 255],   {1:2; 3},    [2; 255])

                % Size-zero bin.
                test(ones(0, 1),   {1:0},         [255]        )
                test([2;3; 255],   {1:2; 1:0; 3}, [3; 255; 255])
            end
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
