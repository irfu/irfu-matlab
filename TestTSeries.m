classdef TestTSeries < matlab.unittest.TestCase
  % Test cases for TSeries and derived classes
  properties
    data1;
    data2;
    time1;
    time2;
    scalarTS1;
    scalarTS2;
  end
  methods (TestClassSetup)
    function createTS(testCase)
      testCase.data1 = (0:1:5)'; % Simple seq. of test data
      % A second seq. of test data (distinguishable from data1), negative values
      testCase.data2 = -10 * testCase.data1;
      % Create simple test time (matching length of data1 and data2).
      testCase.time1 = EpochTT('2018-10-25T10:00:00Z'):1:EpochTT('2018-10-25T10:00:05Z');
      testCase.time2 = testCase.time1 + 0.5; % A second test time (0.5 seconds offset from time1)
      % Create simple scalar TSeries
      testCase.scalarTS1 = irf.ts_scalar(testCase.time1, testCase.data1);
      testCase.scalarTS2 = irf.ts_scalar(testCase.time2, testCase.data2);
    end
  end
  methods(Test)
    function testCombineScalarTS(testCase)
      % combine the two TSeries and sort based on time
      % Expect the two times to be interwoven
      expectedTime = EpochTT('2018-10-25T10:00:00Z'):0.5:EpochTT('2018-10-25T10:00:05.5Z');
      % and the data to be interwoven as well.
      expectedData = zeros(length(testCase.data1)+length(testCase.data2),1);
      expectedData(1:2:end-1) = testCase.data1;
      expectedData(2:2:end)   = testCase.data2;
      % scalarTS1 combined with scalarTS2
      resultTS = combine(testCase.scalarTS1, testCase.scalarTS2);
      testCase.verifyEqual(resultTS.time, expectedTime);
      testCase.verifyEqual(resultTS.data, expectedData);
      % scalarTS2 combined with scalarTS1
      resultTS = combine(testCase.scalarTS2, testCase.scalarTS1);
      testCase.verifyEqual(resultTS.time, expectedTime);
      testCase.verifyEqual(resultTS.data, expectedData);
    end

    function testAbsoluteTS(testCase)
      % test absolute value
      % Positive data, should remain positve and equal.
      expectedData = abs(testCase.data1);
      resultTS = abs(testCase.scalarTS1);
      testCase.verifyEqual(resultTS.data, expectedData);
      % Negative data, should simply change sign.
      expectedData = -testCase.data2;
      resultTS = abs(testCase.scalarTS2);
      testCase.verifyEqual(resultTS.data, expectedData);
      % Absolute function should not alter the time (at all)
      expectedTime = testCase.time2;
      testCase.verifyEqual(resultTS.time, expectedTime);
    end

  end
end