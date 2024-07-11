%
% matlab.unittest automatic test code for both
% solo.adm.parse_dataset_filename(), AND
% solo.adm.create_dataset_filename().
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

      function test_parse(filename, expR)
        actR = solo.adm.parse_dataset_filename(filename);
        testCase.assertEqual(actR, expR)
      end
      function test_create(R, expFilename)
        actFilename = solo.adm.create_dataset_filename(R);
        testCase.assertEqual(actFilename, expFilename)
      end

      function test_both(filename, R)
        % NOTE: Tests pair of conversions: one in each direction.
        test_parse(filename, R)
        test_create(R, filename)
      end



      %==========================
      % NOT dataset, but similar
      %==========================
      test_parse('SOLO_HK_RPW-BIA_V01.skt', ...
        []);
      test_parse('SOLO_HK_RPW-BIA_V01.CDF', ...
        []);   % Uppercase suffix.
      test_parse('solo_l2_rpw-lfr-surv-swf-e-cdag_20200213_V01.cdf', ...
        []);   % Lowercase L2.
      test_parse('solo_L3_rpw-lfr-surv-cwf-e_20200423_V01.png', ...
        []);   % Summary plot.

      % Master CDFs
      test_parse('SOLO_HK_RPW-BIA_V01.cdf', ...
        []);   % No time interval. (Master CDF?)
      test_parse('SOLO_L1_RPW-BIA-CURRENT_V01.Draft.cdf', ...
        []);   % No time interval. (Master CDF?)



      %=======================
      % Ground testing, "LES"
      %=======================
      test_both(...
        'solo_HK_rpw-bia_20170531T173115-20170531T173215_V01_les-ec7fe58.cdf', ...
        struct(...
        'datasetId',       'SOLO_HK_RPW-BIA', ...
        'isCdag',          false, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_HK_rpw-bia', ...
        'versionStr',      '01', ...
        'dateVec1',        [2017 05 31 17 31 15], ...
        'dateVec2',        [2017 05 31 17 32 15], ...
        'timeIntervalStr', '20170531T173115-20170531T173215', ...
        'lesTestStr',      'les-ec7fe58', ...
        'unoffExtension',  []));
      test_both(...
        'solo_HK_rpw-bia_20190523T080316-20190523T134337_V01_les-7ae6b5e.ORIGINAL.cdf', ...
        struct(...
        'datasetId',       'SOLO_HK_RPW-BIA', ...
        'isCdag',          false, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_HK_rpw-bia', ...
        'versionStr',      '01', ...
        'dateVec1',        [2019 05 23 08 03 16], ...
        'dateVec2',        [2019 05 23 13 43 37], ...
        'timeIntervalStr', '20190523T080316-20190523T134337', ...
        'lesTestStr',      'les-7ae6b5e', ...
        'unoffExtension',  {{'ORIGINAL'}}));

      %=======================
      % Ground testing, "CNE"
      %=======================
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.cdf', ...
        struct(...
        'datasetId',       'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',          false, ...
        'dsicdagCase',     'upper', ...
        'fnDatasetIdCdag', 'SOLO_L1_RPW-BIA-CURRENT', ...
        'versionStr',      '02', ...
        'cneTestStr',      '4129f0b_CNE', ...
        'unoffExtension', []));
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.modified.cdf', ...
        struct(...
        'datasetId',       'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',          false, ...
        'dsicdagCase',     'upper', ...
        'fnDatasetIdCdag', 'SOLO_L1_RPW-BIA-CURRENT', ...
        'versionStr',      '02', ...
        'cneTestStr',      '4129f0b_CNE', ...
        'unoffExtension',  {{'modified'}}));

      %=================
      % In-flight, date
      %=================
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e_20200213_V01.cdf', ...
        struct(...
        'datasetId',       'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',          false, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_L2_rpw-lfr-surv-swf-e', ...
        'versionStr',      '01', ...
        'dateVec',         [2020 02 13], ...
        'timeIntervalStr', '20200213', ...
        'unoffExtension',  []));
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e-cdag_20200213_V012345.MODIF.cdf', ...
        struct(...
        'datasetId',       'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',          true, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_L2_rpw-lfr-surv-swf-e-cdag', ...
        'versionStr',      '012345', ...
        'dateVec',         [2020 02 13], ...
        'timeIntervalStr', '20200213', ...
        'unoffExtension',  {{'MODIF'}}));    % Test longer version string.

      % NOTE: MAG
      test_both(...
        'solo_L2_mag-srf-normal_20201231_V03.cdf', ...
        struct(...
        'datasetId',       'SOLO_L2_MAG-SRF-NORMAL', ...
        'isCdag',          false, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_L2_mag-srf-normal', ...
        'versionStr',      '03', ...
        'dateVec',         [2020 12 31], ...
        'timeIntervalStr', '20201231', ...
        'unoffExtension',  []));

      % NOTE: RCT / level "CAL"
      test_both(...
        'solo_CAL_rpw-bias_20220210-20990101_V01.cdf', ...
        struct(...
        'datasetId',       'SOLO_CAL_RPW-BIAS', ...
        'isCdag',          false, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_CAL_rpw-bias', ...
        'versionStr',      '01', ...
        'dateVec1',         [2022 02 10], ...
        'dateVec2',         [2099 01 01], ...
        'timeIntervalStr', '20220210-20990101', ...
        'unoffExtension',  []));    % Test longer version string.

      %==========================
      % In-flight, time interval
      %==========================
      test_both(...
        'solo_L1_rpw-bia-sweep-cdag_20190102T030405-20201112T131415_V01.cdf', ...
        struct(...
        'datasetId',       'SOLO_L1_RPW-BIA-SWEEP', ...
        'isCdag',          true, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_L1_rpw-bia-sweep-cdag', ...
        'versionStr',      '01', ...
        'dateVec1',        [2019 01 02 03 04 05], ...
        'dateVec2',        [2020 11 12 13 14 15], ...
        'timeIntervalStr', '20190102T030405-20201112T131415', ...
        'unoffExtension',  []));   % Guiness world record-length sweep...
      test_both(...
        'solo_L1_rpw-bia-current-cdag_20200401-20200421_V01.cdf', ...
        struct(...
        'datasetId',       'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',          true, ...
        'dsicdagCase',     'lower', ...
        'fnDatasetIdCdag', 'solo_L1_rpw-bia-current-cdag', ...
        'versionStr',      '01', ...
        'dateVec1',        [2020 04 01], ...
        'dateVec2',        [2020 04 21], ...
        'timeIntervalStr', '20200401-20200421', ...
        'unoffExtension',  []));

    end



  end    % methods(Test)

end
