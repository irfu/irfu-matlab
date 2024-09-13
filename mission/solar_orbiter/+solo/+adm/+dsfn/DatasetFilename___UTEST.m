%
% matlab.unittest automatic test code for solo.adm.dsfn.DatasetFilename.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DatasetFilename___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Test more filename non-matches.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test_parse_filename(filename, ExpDfStruct)
        ActDf = solo.adm.dsfn.DatasetFilename.parse_filename(filename);

        if isstruct(ExpDfStruct)
          % Remove private property.
          ActDfStruct = struct(ActDf);
          ActDfStruct = rmfield(ActDfStruct, 'dsicdagUppercase');

          ExpDfStruct.filename = filename;

          % NOTE: struct(object) also returns fields for private properties!
          testCase.assertEqual(ActDfStruct, ExpDfStruct)

        elseif isempty(ExpDfStruct)
          testCase.assertEqual(ActDf, ExpDfStruct)

        end
      end

      function test_filename(S, expFilename)
        Df = solo.adm.dsfn.DatasetFilename(S);
        actFilename = Df.filename;

        testCase.assertEqual(actFilename, expFilename)
      end

      function test_both(filename, DfStruct)
        % NOTE: Tests pair of conversions: one in each direction.
        test_parse_filename(filename, DfStruct)
        test_filename(DfStruct, filename)
      end

      NAT = irf.dt.UTC('NaT');



      %=====================================================
      % NOT dataset filenames, but similar (all unsupported)
      %=====================================================
      test_parse_filename('SOLO_HK_RPW-BIA_V01.skt', ...
        []);
      test_parse_filename('SOLO_HK_RPW-BIA_V01.CDF', ...
        []);   % Uppercase suffix, no time interval string.
      test_parse_filename('solo_l2_rpw-lfr-surv-swf-e-cdag_20200213_V01.cdf', ...
        []);   % Lowercase L2.
      test_parse_filename('solo_L3_rpw-lfr-surv-cwf-e_20200423_V01.png', ...
        []);   % Summary plot.

      % Master CDFs
      test_parse_filename('SOLO_HK_RPW-BIA_V01.cdf', ...
        []);   % No time interval. (Master CDF?)
      test_parse_filename('SOLO_L1_RPW-BIA-CURRENT_V01.Draft.cdf', ...
        []);   % No time interval. (Master CDF?)



      %=======================
      % Ground testing, "LES"
      %=======================
      test_both(...
        'solo_HK_rpw-bia_20170531T173115-20170531T173215_V01_les-ec7fe58.cdf', ...
        struct(...
        'datasetId',          'SOLO_HK_RPW-BIA', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'solo_HK_rpw-bia', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2017 05 31 17 31 15]), ...
        'Dt2',                irf.dt.UTC([2017 05 31 17 32 15]), ...
        'timeIntervalFormat', 'SECOND_TO_SECOND', ...
        'timeIntervalStr',    '20170531T173115-20170531T173215', ...
        'lesTestStr',         'les-ec7fe58', ...
        'cneTestStr',         [] ...
        ))
      test_both(...
        'solo_HK_rpw-bia_20190523T080316-20190523T134337_V01_les-7ae6b5e.cdf', ...
        struct(...
        'datasetId',          'SOLO_HK_RPW-BIA', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'solo_HK_rpw-bia', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2019 05 23 08 03 16]), ...
        'Dt2',                irf.dt.UTC([2019 05 23 13 43 37]), ...
        'timeIntervalFormat', 'SECOND_TO_SECOND', ...
        'timeIntervalStr',    '20190523T080316-20190523T134337', ...
        'lesTestStr',         'les-7ae6b5e', ...
        'cneTestStr',         []))

      %=======================
      % Ground testing, "CNE"
      %=======================
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'SOLO_L1_RPW-BIA-CURRENT', ...
        'versionNbr',         2, ...
        'versionStr',         '02', ...
        'Dt1',                NAT, ...
        'Dt2',                NAT, ...
        'timeIntervalFormat', 'NO_TIME_INTERVAL', ...
        'timeIntervalStr',    [], ...
        'lesTestStr',         [], ...
        'cneTestStr',         '4129f0b_CNE'))
      test_both(...
        'SOLO_L1_RPW-BIA-CURRENT_4129f0b_CNE_V02.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'SOLO_L1_RPW-BIA-CURRENT', ...
        'versionNbr',         2, ...
        'versionStr',        '02', ...
        'Dt1',                NAT, ...
        'Dt2',                NAT, ...
        'timeIntervalFormat', 'NO_TIME_INTERVAL', ...
        'timeIntervalStr',    [], ...
        'lesTestStr',         [], ...
        'cneTestStr',         '4129f0b_CNE'))


      %===========
      % In-flight
      %===========
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e_20200213_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'solo_L2_rpw-lfr-surv-swf-e', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2020 02 13, 00 00 00]), ...
        'Dt2',                irf.dt.UTC([2020 02 14, 00 00 00]), ...
        'timeIntervalFormat', 'DAY', ...
        'timeIntervalStr',    '20200213', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))

      % Test longer version string.
      test_both(...
        'solo_L2_rpw-lfr-surv-swf-e-cdag_20200213_V12345.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
        'isCdag',             true, ...
        'filenameDsiCdag',    'solo_L2_rpw-lfr-surv-swf-e-cdag', ...
        'versionNbr',         12345, ...
        'versionStr',         '12345', ...
        'Dt1',                irf.dt.UTC([2020 02 13 00 00 00]), ...
        'Dt2',                irf.dt.UTC([2020 02 14 00 00 00]), ...
        'timeIntervalFormat', 'DAY', ...
        'timeIntervalStr',    '20200213', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))

      % NOTE: MAG
      test_both(...
        'solo_L2_mag-srf-normal_20201231_V03.cdf', ...
        struct(...
        'datasetId',          'SOLO_L2_MAG-SRF-NORMAL', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'solo_L2_mag-srf-normal', ...
        'versionNbr',         3, ...
        'versionStr',         '03', ...
        'Dt1',                irf.dt.UTC([2020 12 31 00 00 00]), ...
        'Dt2',                irf.dt.UTC([2021 01 01 00 00 00]), ...
        'timeIntervalFormat', 'DAY', ...
        'timeIntervalStr',    '20201231', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))



      % NOTE: RCT / level "CAL"
      test_both(...
        'solo_CAL_rpw-bias_20220210-20990101_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_CAL_RPW-BIAS', ...
        'isCdag',             false, ...
        'filenameDsiCdag',    'solo_CAL_rpw-bias', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2022 02 10 00 00 00]), ...
        'Dt2',                irf.dt.UTC([2099 01 02 00 00 00]), ...
        'timeIntervalFormat', 'DAY_TO_DAY', ...
        'timeIntervalStr',    '20220210-20990101', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))

      test_both(...
        'solo_L1_rpw-bia-sweep-cdag_20190102T030405-20201112T131415_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-SWEEP', ...
        'isCdag',             true, ...
        'filenameDsiCdag',    'solo_L1_rpw-bia-sweep-cdag', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2019 01 02 03 04 05]), ...
        'Dt2',                irf.dt.UTC([2020 11 12 13 14 15]), ...
        'timeIntervalFormat', 'SECOND_TO_SECOND', ...
        'timeIntervalStr',    '20190102T030405-20201112T131415', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))

      test_both(...
        'solo_L1_rpw-bia-current-cdag_20200401-20200421_V01.cdf', ...
        struct(...
        'datasetId',          'SOLO_L1_RPW-BIA-CURRENT', ...
        'isCdag',             true, ...
        'filenameDsiCdag',    'solo_L1_rpw-bia-current-cdag', ...
        'versionNbr',         1, ...
        'versionStr',         '01', ...
        'Dt1',                irf.dt.UTC([2020 04 01 00 00 00]), ...
        'Dt2',                irf.dt.UTC([2020 04 22 00 00 00]), ...
        'timeIntervalFormat', 'DAY_TO_DAY', ...
        'timeIntervalStr',    '20200401-20200421', ...
        'lesTestStr',         [], ...
        'cneTestStr',         []))
    end



  end    % methods(Test)



end