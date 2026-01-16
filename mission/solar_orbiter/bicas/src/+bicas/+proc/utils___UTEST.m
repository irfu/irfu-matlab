%
% matlab.unittest automatic test code for bicas.proc.utils().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08 from older test code.
%
classdef utils___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_assert_increasing(T)
      function test(array, isMonotonic)
        errorId = 'test:Error';
        msg = '<Error message>';

        % Test both with and without transposition.
        bicas.proc.utils.assert_increasing(array, isMonotonic, errorId, msg);
        bicas.proc.utils.assert_increasing(array', isMonotonic, errorId, msg);
      end

      function test_exc(array, isMonotonic)
        errorId = 'test:Error';
        msg = '<Error message>';
        T.verifyError(...
          @() bicas.proc.utils.assert_increasing(array, isMonotonic, errorId, msg), ...
          ?MException)
      end

      for isMonotonic = [false, true]
        test(zeros(0, 1), isMonotonic)
        test([3], isMonotonic)
        test([1,2,3,4,5], isMonotonic)

        test_exc([5,4,6,7,8], isMonotonic)
      end

      test([1,2,3,3,4], false)
      test_exc([1,2,3,3,4], true)
    end



    function test_set_NaN_end_of_rows(T)

      function test(zv, snapshotLengths, expZv)
        actZv = bicas.proc.utils.set_NaN_end_of_rows(zv, snapshotLengths);
        T.verifyEqual(actZv, expZv)
      end

      test(ones(0,4),              ones(0,1), ones(0,4));
      test([0,1,2],                [3],       [0,1,2]);
      test([0,1,2,3,4; 5,6,7,8,9], [2;4],     [0,1,NaN,NaN,NaN; 5,6,7,8,NaN]);
    end



    function test_find_data_gaps(T)
      function test(tt2000Ar, samplRateHz, maxSampleGapFraction, expISegmentAr)
        tt2000Ar = int64(tt2000Ar);

        actISegmentAr = bicas.proc.utils.find_data_gaps(...
          tt2000Ar, samplRateHz, maxSampleGapFraction);

        T.assertEqual(actISegmentAr, expISegmentAr)
      end

      % Zero data points
      test(zeros(0, 1), zeros(0, 1), 2.01, zeros(0, 1))

      % One data point
      test([10] * 1e9, [3  ], 2.01, [0])
      test([10] * 1e9, [0.1], 0.01, [0])

      % Several data points
      % One sampling rate.
      % No data gap.
      test(...
        [0   10   20   30  ]' * 1e9, ...
        [0.1  0.1  0.1  0.1]', 1.01, ...
        [0    0    0    0  ]')
      test(...
        [0   10   20   30  ]' * 1e9, ...
        [0.1  0.1  0.1  0.1]', 0.99, ...
        [0    1    2    3  ]')

      % One sampling rate, one data gap.
      test(...
        [0 1 2   4 5 6 7]' * 1e9, ...
        [1 1 1   1 1 1 1]', 2.01, ...
        [0 0 0   0 0 0 0]')
      test(...
        [0 1 2   4 5 6 7]' * 1e9, ...
        [1 1 1   1 1 1 1]', 1.99, ...
        [0 0 0   1 1 1 1]')

      % Two sampling rates.
      % No data gap.
      test(...
        [2   4     7    11   ]' * 1e9, ...
        [0.5 0.5   0.25  0.25]', 1.01, ...
        [0   0     0     0   ]')
      test(...
        [2   4     7    11   ]' * 1e9, ...
        [0.5 0.5   0.25  0.25]', 0.99, ...
        [0   1     2     3   ]')

      % Two sampling rates
      % One data gap (small).
      test(...
        [0 1 2   4   6   8  ]' * 1e9, ...
        [1 1 1   0.5 0.5 0.5]', 2/(0.5+1) * 1.01, ...
        [0 0 0   0   0   0]')
      test(...
        [0 1 2   4   6   8  ]' * 1e9, ...
        [1 1 1   0.5 0.5 0.5]', 2/(0.5+1) * 0.99, ...
        [0 0 0   1   1   1]')

      % Complex, "realistic" test
      % Fluctuating sampling rate in timestamps.
      test(...
        [0.0 2.1 4.0 6.0  10.0 11.1 12.0 13.1   20   30   40  ]' * 1e9, ...
        [0.5 0.5 0.5 0.5   1    1    1    1      0.1  0.1  0.1]', 1.1, ...
        [0   0   0   0     1    1    1    1      2    2    2  ]')
    end



    function test_convert_matrix_to_cell_array_of_vectors(T)

      function test(M, nCopyColsPerRowArray, expCa)
        actCa = bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
          M, nCopyColsPerRowArray);

        T.verifyEqual(actCa, expCa)
      end

      test(zeros(0,1),  zeros(0,1), cell(0,1));

      test([1,2,3,4,5], [0],        {zeros(0, 1)});
      test([1,2,3,4,5], [3],        {[1 2 3]'});
      test([1,2,3,4,5], [5],        {[1 2 3 4 5]'});

      test([1,2,3,4,5; 6,7,8,9,0], [3; 2], {[1 2 3]'; [6 7]'});
      test([1,2,3,4,5; 6,7,8,9,0], [0; 5], {zeros(0, 1); [6 7 8 9 0]'});
    end



    function test_convert_cell_array_of_vectors_to_matrix(T)

      function test(ca, nMatrixColumns, expM, expNCopyColsPerRowVec)

        [actM, actNCopyColsPerRowVec] = ...
          bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
          ca, nMatrixColumns);
        T.verifyEqual(expM,                  actM)
        T.verifyEqual(expNCopyColsPerRowVec, actNCopyColsPerRowVec)
      end

      % Zero rows
      for nCols = [0, 10]
        test(...
          cell(0, 1), nCols, ...
          zeros(0, nCols), zeros(0, 1) ...
          )
      end

      % One row
      test(...
        {[]}, 5, ...
        [nan, nan, nan, nan, nan], [0] ...
        )
      test(...
        {[1, 2, 3]}, 5, ...
        [1, 2, 3, nan, nan], [3] ...
        )

      % >1 rows
      test(...
        {[1, 2, 3]; [11, 12]; []}, 3, ...
        [1,2,3; 11,12,nan; nan,nan,nan], [3; 2; 0] ...
        )
    end



  end    % methods(Test)



end
