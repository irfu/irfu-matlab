%
% matlab.unittest automatic test code for both
% solo.adm.dsfn.parse_dataset_filename(), AND
% solo.adm.dsfn.create_dataset_filename().
%
%
% NOTE: Tests TWO functions which are each other's inverses.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef parse_dataset_filename___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      % PROPOSAL: Test more non-matches.

      function test_parse(filename, expR, expS)
        [actR, actS] = solo.adm.dsfn.parse_dataset_filename(filename);

        testCase.assertEqual(actR, expR)
        testCase.assertEqual(actS, expS)
      end

      function test_create(R, expFilename)
        actFilename = solo.adm.dsfn.create_dataset_filename(R);
        testCase.assertEqual(actFilename, expFilename)
      end

      function test_both(filename, R, expS)
        % NOTE: Tests pair of conversions: one in each direction.
        test_parse(filename, R, expS)
        test_create(R, filename)
      end



      %=====================================================
      % NOT dataset filenames, but similar (all unsupported)
      %=====================================================
      test_parse('SOLO_HK_RPW-BIA_V01.skt', ...
        [], []);
      test_parse('SOLO_HK_RPW-BIA_V01.CDF', ...
        [], []);   % Uppercase suffix, no time interval string.
      test_parse('solo_l2_rpw-lfr-surv-swf-e-cdag_20200213_V01.cdf', ...
        [], []);   % Lowercase L2.
      test_parse('solo_L3_rpw-lfr-surv-cwf-e_20200423_V01.png', ...
        [], []);   % Summary plot.

      % Master CDFs
      test_parse('SOLO_HK_RPW-BIA_V01.cdf', ...
        [], []);   % No time interval. (Master CDF?)
      test_parse('SOLO_L1_RPW-BIA-CURRENT_V01.Draft.cdf', ...
        [], []);   % No time interval. (Master CDF?)



      %=======================
      % Ground testing, "LES"
      %=======================
      test_both(...
        'solo_HK_rpw-bia_20170531T173115-20170531T173215_V01_les-ec7fe58.cdf', ...
        struct(...
        'datasetId',          'SOLO_HK_RPW-BIA', ...
        'isCdag',             false, ...
        'versionStr',         '01', ...
        'dateVec1',           [2017 05 31 17 31 15], ...
        'dateVec2',           [2017 05 31 17 32 15], ...
        'timeIntervalFormat', 'SECOND_TO_SECOND', ...
        'lesTestStr',         'les-ec7fe58'), ...
        struct(...
        'timeIntervalStr', '20170531T173115-20170531T173215', ...
        'filenameDsiCdag', 'solo_HK_rpw-bia'));
      test_both(...
        'solo_HK_rpw-bia_20190523T080316-20190523T134337_V01_les-7ae6b5e.cdf', ...
        struct(...
        'datasetId',          'SOLO_HK_RPW-BIA', ...
        'isCdag',             false, ...
        'versionStr',         '01', ...
        'dateVec1',           [2019 05 23 08 03 16], ...
        'dateVec2',           [2019 05 23 13 43 37], ...
        'timeIntervalFormat', 'SECOND_TO_SECOND', ...
        'lesTestStr',         'les-7ae6b5e'), ...
        struct(...
        'timeIntervalStr', '20190523T080316-20190523T134337', ...
        'filenameDsiCdag', 'solo_HK_rpw-bia'));

      %=======================
      % Ground testing, "CNE"
      %=======================
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             false, ...
        'versionStr',         '02', ...
        'dateVec1',           [0 0 0 0 0 0], ...
        'dateVec2',           [0 0 0 0 0 0], ...
        'timeIntervalFormat', 'NO_TIME_INTERVAL', ...
        'cneTestStr',         '4129f0b_CNE'), ...
        struct(...
        'timeIntervalStr', [], ...
        'filenameDsiCdag', 'SOLO_L1_RPW-BIA-CURRENT'));
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             false, ...
        'versionStr',         '02', ...
        'dateVec1',           [0 0 0 0 0 0], ...
        'dateVec2',           [0 0 0 0 0 0], ...
        'timeIntervalFormat', 'NO_TIME_INTERVAL', ...
        'cneTestStr',         '4129f0b_CNE'), ...
        struct(...
        'timeIntervalStr', [], ...
        'filenameDsiCdag', 'SOLO_L1_RPW-BIA-CURRENT'));

      %===========
      % In-flight
      %===========
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e_20200213_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',             false, ...
        'versionStr',         '01', ...
        'dateVec1',           [2020 02 13, 00 00 00], ...
        'dateVec2',           [2020 02 14, 00 00 00], ...
        'timeIntervalFormat', 'DAY'), ...
        struct(...
        'timeIntervalStr', '20200213', ...
        'filenameDsiCdag', 'solo_L2_rpw-lfr-surv-swf-e'));

      % Test longer version string.
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e-cdag_20200213_V012345.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',             true, ...
        'versionStr',         '012345', ...
        'dateVec1',           [2020 02 13 00 00 00], ...
        'dateVec2',           [2020 02 14 00 00 00], ...
        'timeIntervalFormat', 'DAY'), ...
        struct(...
        'timeIntervalStr', '20200213', ...
        'filenameDsiCdag', 'solo_L2_rpw-lfr-surv-swf-e-cdag'))

      % NOTE: MAG
      test_both(...
        'solo_L2_mag-srf-normal_20201231_V03.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_MAG-SRF-NORMAL', ...
        'isCdag',             false, ...
        'versionStr',         '03', ...
        'dateVec1',           [2020 12 31 00 00 00], ...
        'dateVec2',           [2021 01 01 00 00 00], ...
        'timeIntervalFormat', 'DAY'), ...
        struct(...
        'timeIntervalStr', '20201231', ...
        'filenameDsiCdag', 'solo_L2_mag-srf-normal'))

      % NOTE: RCT / level "CAL"
      test_both(...
        'solo_CAL_rpw-bias_20220210-20990101_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_CAL_RPW-BIAS', ...
        'isCdag',             false, ...
        'versionStr',         '01', ...
        'dateVec1',           [2022 02 10 00 00 00], ...
        'dateVec2',           [2099 01 02 00 00 00], ...
        'timeIntervalFormat', 'DAY_TO_DAY'), ...
        struct(...
        'timeIntervalStr', '20220210-20990101', ...
        'filenameDsiCdag', 'solo_CAL_rpw-bias'));

      test_both(...
        'solo_L1_rpw-bia-sweep-cdag_20190102T030405-20201112T131415_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-SWEEP', ...
        'isCdag',             true, ...
        'versionStr',         '01', ...
        'dateVec1',           [2019 01 02 03 04 05], ...
        'dateVec2',           [2020 11 12 13 14 15], ...
        'timeIntervalFormat', 'SECOND_TO_SECOND'), ...
        struct(...
        'timeIntervalStr', '20190102T030405-20201112T131415', ...
        'filenameDsiCdag', 'solo_L1_rpw-bia-sweep-cdag'))

      test_both(...
        'solo_L1_rpw-bia-current-cdag_20200401-20200421_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             true, ...
        'versionStr',         '01', ...
        'dateVec1',           [2020 04 01 00 00 00], ...
        'dateVec2',           [2020 04 22 00 00 00], ...
        'timeIntervalFormat', 'DAY_TO_DAY'), ...
        struct(...
        'timeIntervalStr', '20200401-20200421', ...
        'filenameDsiCdag', 'solo_L1_rpw-bia-current-cdag'))
    end



  end    % methods(Test)

end
