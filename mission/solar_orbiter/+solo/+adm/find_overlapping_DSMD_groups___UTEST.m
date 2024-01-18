%
% matlab.unittest automatic test code for
% solo.adm.find_overlapping_DSMD_groups().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef find_overlapping_DSMD_groups___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      % PROPOSAL: Abolish function Dsmd(). Convert filenames-->DSMD in
      %           test().

      function test(dsmdArray, datasetIdList, expDsmdGroupsCa)
        % dsmdArray       = solo.adm.paths_to_DSMD_array(filenamesCa);
        % expDsmdGroupsCa = solo.adm.paths_to_DSMD_array(expFilenamesCa);

        actDsmdGroupsCa = solo.adm.find_overlapping_DSMD_groups(...
          dsmdArray, datasetIdList);

        testCase.assertEqual(actDsmdGroupsCa, expDsmdGroupsCa)
        testCase.assertTrue(iscell(expDsmdGroupsCa))
        testCase.assertTrue(iscolumn(expDsmdGroupsCa))
      end

      function Dsmd = DSMD(filename)
        irf.assert.castring(filename);
        Dsmd = solo.adm.paths_to_DSMD_array({filename});
      end

      % SW_13    = DSMD('solo_L1_rpw-bia-sweep-cdag_20200313T100000-20200313T110000_V01.cdf');
      % SW_14    = DSMD('solo_L1_rpw-bia-sweep-cdag_20200314T100000-20200314T110000_V01.cdf');
      SW_13_14 = DSMD('solo_L1_rpw-bia-sweep-cdag_20200313T100000-20200314T110000_V01.cdf');

      HK_13    = DSMD('solo_HK_rpw-bia_20200313_V05.cdf');
      HK_14    = DSMD('solo_HK_rpw-bia_20200314_V05.cdf');
      HK_15    = DSMD('solo_HK_rpw-bia_20200315_V05.cdf');

      % CWF_13 = DSMD('solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200313_V03.cdf');
      % SWF_13 = DSMD('solo_L1R_rpw-lfr-surv-swf-e-cdag_20200313_V03.cdf');
      CWF_14 = DSMD('solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200314_V03.cdf');
      SWF_14 = DSMD('solo_L1R_rpw-lfr-surv-swf-e-cdag_20200314_V03.cdf');
      % CWF_15 = DSMD('solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200315_V03.cdf');
      % SWF_15 = DSMD('solo_L1R_rpw-lfr-surv-swf-e-cdag_20200315_V03.cdf');

      % CUR_02 = DSMD('solo_L1_rpw-bia-current-cdag_20200201T000000-20200301T000000_V01.cdf');
      % CUR_03 = DSMD('solo_L1_rpw-bia-current-cdag_20200301T000000-20200401T000000_V01.cdf');

      % NOTE: Implementation has special case for empty array.
      test(solo.adm.DSMD.empty(1, 0), {'SOLO_HK_RPW-BIA'}, cell(0, 1));

      test([HK_14], {'SOLO_HK_RPW-BIA'}, {HK_14});
      test([HK_13, HK_15], {'SOLO_HK_RPW-BIA'}, {HK_13; HK_15});
      test([HK_13, HK_14], {'SOLO_HK_RPW-BIA'}, {HK_13; HK_14});

      test([HK_14, CWF_14], {'SOLO_HK_RPW-BIA', 'SOLO_L1R_RPW-LFR-SURV-CWF-E'}, {[HK_14; CWF_14]});
      test([HK_14, CWF_14], {'SOLO_L1R_RPW-LFR-SURV-CWF-E', 'SOLO_HK_RPW-BIA'}, {[HK_14; CWF_14]});
      test([HK_14, CWF_14, SWF_14], ...
        {'SOLO_HK_RPW-BIA', 'SOLO_L1R_RPW-LFR-SURV-CWF-E', 'SOLO_L1R_RPW-LFR-SURV-SWF-E'}, ...
        {[HK_14; CWF_14; SWF_14]});

      test([HK_13, HK_14, HK_15, CWF_14], {'SOLO_HK_RPW-BIA', 'SOLO_L1R_RPW-LFR-SURV-CWF-E'}, ...
        {[HK_14; CWF_14]});
      test([HK_14, CWF_14, SW_13_14], {'SOLO_HK_RPW-BIA', 'SOLO_L1R_RPW-LFR-SURV-CWF-E'}, ...
        {[HK_14; CWF_14]});
      test([HK_13, HK_14,  SW_13_14, HK_15], {'SOLO_HK_RPW-BIA', 'SOLO_L1_RPW-BIA-SWEEP'}, ...
        {[HK_13; SW_13_14]; [HK_14; SW_13_14]});
    end



  end    % methods(Test)

end
