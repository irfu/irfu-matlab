%
% Automatic test code.
%
% NOTE: Very low code coverage.
%
function proc_sub23___ATEST()
    downsample_sci_zVar___ATEST()
    %downsample_bin_sci_values___ATEST()
    downsample_Epoch___ATEST
end



% Internal utility function
function assert_iRecordsDwnsCa(iRecordsDwnsCa, zv)
    
    % ~Normalize
    iRecordsDwnsNormCa = arrayfun(...
        @(s) (s{1}(:)'), ...
        iRecordsDwnsCa, ...
        'UniformOutput', false);
    
    iRecordsArray = [iRecordsDwnsNormCa{:}];
    
    assert(all(unique(iRecordsArray) == 1:size(zv, 1)))
end



function downsample_sci_zVar___ATEST
    
    tl = {};
    
    function add_test_1_bin(zv, nMinReqSamples, med, mstd)
        tl{end+1} = EJ_library.atest.CompareFuncResult(...
            @bicas.proc_sub23.downsample_sci_zVar, ...
            {zv, nMinReqSamples, {1:size(zv,1)}}, ...
            {med, mstd});
    end

    function add_test_N_bin(zv, nMinReqSamples, iRecordsDwnsCa, med, mstd)
        assert(isrow(iRecordsDwnsCa))
        assert_iRecordsDwnsCa(iRecordsDwnsCa, zv)
        
        tl{end+1} = EJ_library.atest.CompareFuncResult(...
            @bicas.proc_sub23.downsample_sci_zVar, ...
            {zv, nMinReqSamples, iRecordsDwnsCa'}, ...
            {med, mstd});
    end

    % ERA = Empty Rows (size=0) Array
    ERA = zeros(1,0);

    
    if 1
    add_test_N_bin([1,2; 2,3], 1, {1, 2}, [1,2; 2,3], NaN(2,2))
    
    add_test_N_bin(...
        [1, 2; ...
         2, 3; ...
         3, 4; ...
         4, 5; ...
         5, 6], 1, {1:2, 3:5}, ...
        [1.5, 2.5; ...
         4,   5   ], ...
        [sqrt(2*0.5^2)/1 * [1,1]; ...
         sqrt(4*1^2  )/2 * [1,1]])
     
    add_test_N_bin(...
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
    add_test_N_bin([1,2; 2,3; 3,4; 4,5; 5,6], 3, {1:2, [], 3:5}, ...
        [NaN, NaN; ...
         NaN, NaN; ...
         4, 5], ...
        [NaN, NaN; ...
         NaN, NaN; 
         sqrt(4*1^2  )/2 * [1,1]])
    end
    
    
    
    if 1
    % Empty data
    add_test_1_bin(zeros(0,0), 0, ERA, ERA);
    add_test_1_bin(zeros(0,2), 0, [NaN NaN], [NaN, NaN]);
    
    % mstd=0
    add_test_1_bin([1,2,3              ], 0, [1,2,3], [nan,nan,nan]);
    add_test_1_bin([1,2,3; 1,2,3       ], 0, [1,2,3], [0,0,0]);
    add_test_1_bin([1,2,3; 1,2,3; 1,2,3], 0, [1,2,3], [0,0,0]);
    add_test_1_bin([1    ; 1    ; 1    ], 0, [1],     [0]);
    
    % Test nMinReqSamples
    add_test_1_bin([1,2,3; 1,2,3; 1,2,3], 3, [1,2,3], [0,0,0]);
    add_test_1_bin([1,2,3; 1,2,3; 1,2,3], 4, [nan,nan,nan], [nan,nan,nan]);
    
    
    % Average of two values (special case)
    add_test_1_bin([1,2,3; 2,3,4], 0, [1.5, 2.5, 3.5], sqrt(0.5)*[1,1,1]);
    % Nominal median
    add_test_1_bin([1;2;10],       0, [2], sqrt( (1^2+0^2+8^2)/2 ));
    end
    
    
    EJ_library.atest.run_tests(tl)
end


% 
% function downsample_bin_sci_values___ATEST()
%     
%     tl = {};
%     
%     function add_test(zvSegment, nMinReqSamples, med, mstd)
%         tl{end+1} = EJ_library.atest.CompareFuncResult(...
%             @bicas.proc_sub23.downsample_bin_sci_values, ...
%             {zvSegment, nMinReqSamples}, ...
%             {med, mstd});
%     end
%     
% %     newTestExc = @(zVarSegment, nMinReqSamples) (...
% %         EJ_library.atest.CompareFuncResult(...
% %         @bicas.proc_sub23.downsample_bin_sci_values, ...
% %         {zVarSegment}, ...
% %         {med, mstd}));
% 
%     ERA = zeros(1,0);
%     
%     % Empty data
%     add_test(zeros(0,0), 0, ERA, ERA);
%     add_test(zeros(0,2), 0, [NaN NaN], [NaN, NaN]);
%     
%     % mstd=0
%     add_test([1,2,3              ], 0, [1,2,3], [nan,nan,nan]);
%     add_test([1,2,3; 1,2,3       ], 0, [1,2,3], [0,0,0]);
%     add_test([1,2,3; 1,2,3; 1,2,3], 0, [1,2,3], [0,0,0]);
%     add_test([1    ; 1    ; 1    ], 0, [1],     [0]);
%     
%     % Test nMinReqSamples
%     add_test([1,2,3; 1,2,3; 1,2,3], 3, [1,2,3], [0,0,0]);
%     add_test([1,2,3; 1,2,3; 1,2,3], 4, [nan,nan,nan], [nan,nan,nan]);
%     
%     
%     % Average of two values (special case)
%     add_test([1,2,3; 2,3,4], 0, [1.5, 2.5, 3.5], sqrt(0.5)*[1,1,1]);
%     % Nominal median
%     add_test([1;2;10],       0, [2], sqrt( (1^2+0^2+8^2)/2 ));
%     
%     EJ_library.atest.run_tests(tl)
% end



function downsample_Epoch___ATEST
    
    function add_test(...
            zvAllUtcCa, ...
            boundaryRefUtc, ...
            intervalLengthWolsNs, ...
            timestampPosWolsNs, ...
            ...
            zvIntervalsUtcCa, ...
            iRecordsCa, ...
            nRecordsPerBin, ...
            binSizeArrayNs)

        tl{end+1} = EJ_library.atest.CompareFuncResult(...
            @bicas.proc_sub23.downsample_Epoch, ...
            ...
            {spdfparsett2000(zvAllUtcCa), ...
            spdfparsett2000(boundaryRefUtc), ...
            int64(intervalLengthWolsNs), ...
            int64(timestampPosWolsNs)}, ...
            ...
            {spdfparsett2000(zvIntervalsUtcCa), ...
            iRecordsCa(:), ...
            nRecordsPerBin(:), ...
            binSizeArrayNs(:)});
    end
    
    
    % ECA = Empty Column (size=0) Array
    ECA = zeros(0,1);
    
    tl = {};
    
    % TEST: No timestamps/empty Epoch.
    add_test({}, '2020-01-01T00:00:00', 10e9, 5e9, {}, {}, [], ECA);
    
    % TEST: One timestamps.
    add_test({'2020-01-01T00:00:06'}, '2020-01-01T00:00:00', 10e9, 5e9, {'2020-01-01T00:00:05'}, {1}, [1], [10e9]);
    add_test({'2020-01-01T00:00:06'}, '2020-01-01T00:00:01', 10e9, 3e9, {'2020-01-01T00:00:04'}, {1}, [1], [10e9]);
    
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
        }, {1, [2,3]'}, [1,2], [10e9, 10e9]);
    
    % TEST: Total interval over LEAP SECOND (end of 2016).
    % TEST: Data gap (two intervals without timestamps).
    % TEST: Boundary reference timestamp far from data.
    % TEST: Generally convoluted case...
    % 
    %tl = {};
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
        }, {[1,2]', [3,4,5,6]', ECA, ECA, [7]}, [2,4,0,0,1], [10e9, 11e9, 10e9, 10e9, 10e9]);
    
    EJ_library.atest.run_tests(tl);
end
