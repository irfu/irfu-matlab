%
% matlab.unittest automatic test code for
% solo.adm.filter_DSMD_CURRENT_largest_time_coverage().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef filter_DSMD_CURRENT_largest_time_coverage___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(inputFileList, expFileList)
        inputDsmdArray = solo.adm.paths_to_DSMD_array( inputFileList);
        expDsmdArray   = solo.adm.paths_to_DSMD_array(expFileList);

        actDsmdArray = solo.adm.filter_DSMD_CURRENT_largest_time_coverage(inputDsmdArray);
        testCase.assertEqual(actDsmdArray, expDsmdArray)
      end

      X1  = 'solo_L2_rpw-lfr-surv-swf-e-cdag_20200314_V01.cdf';
      X2  = 'solo_L2_rpw-lfr-surv-swf-e-cdag_20200314_V02.cdf';

      C1  = 'solo_L1_rpw-bia-current-cdag_20200303T000000-20200321T000000_V01.cdf';
      C2  = 'solo_L1_rpw-bia-current-cdag_20200304T000000-20200401T000000_V01.cdf';
      D1  = 'solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf';
      D21 = 'solo_L1_rpw-bia-current-cdag_20200401T000000-20200501T000000_V01.cdf';
      D22 = 'solo_L1_rpw-bia-current-cdag_20200401T000000-20200501T000000_V02.cdf';   % Same time coverage. Different version.
      %B2 = 'solo_L1_rpw-bia-current-cdag_20200304T000000-20200501T000000_V01.cdf';

      test({C1}, {C1})

      test({D1, D21}, {D21});
      test({C2, C1}, {C2});

      test({C1, D1}, {C1, D1});
      test({C1, D21, C2, D1}, {D21, C2});
      test({C1, D21, X1, C2, D1, X2}, {D21, X1, C2, X2});

      test({C1, D21, D22, X1, C2, D1, X2}, {D21, D22, X1, C2, X2});

      test({D21, D22}, {D21, D22});
      test({D21, D22, D1}, {D21, D22});
      test({D1, D21, D22}, {D21, D22});
    end



  end    % methods(Test)


end
