%
% matlab.unittest automatic test code for
% irf.utils.find_interval_overlaps().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef find_interval_overlaps___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      % Arbitrary number output variables.
      function test(F, expSetsCa, expH)

        expSetsCa = irf.utils.find_interval_overlaps___UTEST.normalize_rv(expSetsCa);
        expH1Array = expH(:,1);
        expH2Array = expH(:,2);
        expNArray = irf.utils.find_interval_overlaps___UTEST.get_n_interval_array(expSetsCa);

        [actSetsCa, actNArray, actH1Array, actH2Array] = ...
          irf.utils.find_interval_overlaps(F(:,1), F(:,2));

        testCase.verifyEqual(actSetsCa, expSetsCa)
        testCase.verifyEqual(actH1Array, expH1Array)
        testCase.verifyEqual(actH2Array, expH2Array)
        testCase.verifyEqual(actNArray, expNArray)
      end

      function test_exc(F)
        testCase.verifyError(...
          @() irf.utils.find_interval_overlaps(...
          F(:,1), F(:,2)), ...
          ?MException...
          )
      end


      % TEMP
      % test([1 2; 1 2],      {[1 2]}, [1, 2]);



      EA = zeros(0,2);

      test(EA,         {}, EA);

      % Negative interval length.
      test_exc([2 1]);

      test([1 2],      {1}, [1, 2]);

      % Intervals with no intervals between them.
      test([1 2; 3 4], {1, [], 2}, [1 2; 2 3; 3 4]);

      % Identical intervals.
      test([1 2; 1 2], {[1 2]}, [1 2]);
      test([3 3; 3 3], {[1 2]}, [3 3]);

      % Overlapping intervals.
      test([1 4; 2 7], {1, [1 2], 2}, [1 2; 2 4; 4 7]);

      % Single ZLI
      test([2 2], {1}, [2 2])

      % Touching NZLIs
      test([2 3; 3 5], {1, [1 2], 2}, [2 3; 3 3; 3 5]);

      % ZLI in the middle of other interval.
      test([2 5; 3 3],      {1, [1, 2], 1},  [2 3; 3 3; 3 5]);

      % ZLI simultaneously in the middle of and adjacent to other
      % interval.
      test([2 3; 3 5; 3 3], {1, [1,2,3], 2}, [2 3; 3 3; 3 5]);

      % Touching NZLIs+ZLI (triple point)
      test([2 3; 3 5; 3 3], {1, [1 2 3], 2}, [2 3; 3 3; 3 5]);


      % Reversely ordered intervals.
      test([4 7; 1 2], {2, [], 1}, [1 2; 2 4; 4 7]);
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)

    function nIntArray = get_n_interval_array(setsCa)
      % Derive the nArray return value from other return value.
      nIntArray = ones(0, 1);
      for i = 1:numel(setsCa)
        nIntArray(i) = numel(setsCa{i});
      end
      nIntArray = nIntArray(:);
    end

    function setsCa = normalize_rv(setsCa)
      for i = 1:numel(setsCa)
        setsCa{i} = setsCa{i}(:);
      end
      setsCa = setsCa(:);
    end

  end    % methods(Static, Access=private)



end
