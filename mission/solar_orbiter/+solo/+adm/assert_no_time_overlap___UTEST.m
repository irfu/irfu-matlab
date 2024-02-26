%
% matlab.unittest automatic test code for
% solo.adm.assert_no_time_overlap().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef assert_no_time_overlap___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(filePathCa)
        DsmdArray = solo.adm.paths_to_DSMD_array(filePathCa(:));

        % Test that does not raise exception.
        solo.adm.assert_no_time_overlap(DsmdArray);
      end

      function test_exc(filePathCa)
        DsmdArray = solo.adm.paths_to_DSMD_array(filePathCa(:));
        testCase.verifyError(...
          @() solo.adm.assert_no_time_overlap(DsmdArray), ...
          ?MException)
      end

      test({})
      test({
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        })
      test({
        'solo_HK_rpw-bia_20200313_V05.cdf', ...
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        })

      % Plain overlap.
      test_exc({
        'solo_HK_rpw-bia_20200314_V04.cdf', ...
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        })

      % Time overlap, but for different dataset IDs.
      test({
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        'solo_HK_rpw-bia_20200315_V05.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200314_V03.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200315_V03.cdf', ...
        })

      test_exc({
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        'solo_HK_rpw-bia_20200315_V05.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200314_V03.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200315_V03.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200316_V03.cdf', ...
        'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200315_V04.cdf', ...
        })

      test_exc({
        'solo_HK_rpw-bia_20200314_V05.cdf', ...
        'solo_HK_rpw-bia_20200315_V05.cdf', ...
        'solo_L1_rpw-bia-current-cdag_20200101-20200120_V01.cdf', ...
        'solo_L1_rpw-bia-current-cdag_20200110-20200130_V21.cdf', ...
        })
    end



  end    % methods(Test)

end
