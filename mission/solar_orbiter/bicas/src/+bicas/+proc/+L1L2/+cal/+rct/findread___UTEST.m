%
% matlab.unittest automatic test code for bicas.proc.L1L2.cal.rct.findread.
%
% NOTE: Only tests one function.
% NOTE: Tests for "bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp()"
% have partly been written in order to try out functionality for testing code
% with file operations.
% NOTE: Creates temporary directory and files for every test, separately.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16
%
classdef findread___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Tests for bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log() for BIAS
  %           RCT. Creates BIAS RCT using bicas.tools.rct.create_RCT() as part of
  %           the test.



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and teardown
    % methods which store/read their own data from the testCase object.
    dir
    L
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      Fixture = testCase.applyFixture(...
        matlab.unittest.fixtures.TemporaryFolderFixture);
      % NOTE: The same fixture should always return the same directory.
      testCase.dir = Fixture.Folder;
      testCase.L   = bicas.Logger('NO_STDOUT', false);
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_BRVF_RCT_path(testCase)
      % Create BRVF.
      ExpDtBegin = datetime('2020-02-10T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDtEnd   = datetime('2100-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      expBiasRctFilename = bicas.tools.rct.create_RCT_filename(ExpDtBegin, ExpDtEnd, 1);
      bicas.tools.rct.create_BRVF(testCase.dir, expBiasRctFilename, ExpDtBegin, ExpDtEnd);

      expBiasRctPath = fullfile(testCase.dir, expBiasRctFilename);
      expBrvfPath    = fullfile(testCase.dir, bicas.const.BRVF_FILENAME);

      %========
      % TEST 1
      %========
      [actBiasRctPath, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.get_BRVF_RCT_path(testCase.dir, ExpDtBegin, ExpDtEnd);

      testCase.assertEqual(actBiasRctPath, expBiasRctPath)
      testCase.assertEqual(actBrvfPath,    expBrvfPath)



      %========
      % TEST 2
      %========
      Duration = duration('00:00:10');
      % NOTE: Only checks for non-error time interval case.
      [actBiasRctPath, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.get_BRVF_RCT_path(...
        testCase.dir, ...
        ExpDtBegin + Duration, ...
        ExpDtEnd   - Duration);
    end



    function test_read_BRVF(testCase)
      % Create BRVF.
      ExpDtBegin = datetime('2020-02-10T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDtEnd   = datetime('2100-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      expBiasRctFilename = bicas.tools.rct.create_RCT_filename(ExpDtBegin, ExpDtEnd, 1);
      bicas.tools.rct.create_BRVF(testCase.dir, expBiasRctFilename, ExpDtBegin, ExpDtEnd);

      [actRctFilename, ActDtValidityBegin, ActDtValidityEnd, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.read_BRVF(testCase.dir);

      testCase.assertEqual(actRctFilename,     expBiasRctFilename)
      testCase.assertEqual(ActDtValidityBegin, ExpDtBegin)
      testCase.assertEqual(ActDtValidityEnd,   ExpDtEnd)
    end



    function test_get_RCT_path_by_regexp_empty(testCase)
      bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(testCase, {});

      testCase.verifyError(...
        @() bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        testCase.dir, '20[0-9][0-9]\.cdf', testCase.L), ...
        'BICAS:CannotFindRegexMatchingRCT')

    end



    function test_get_RCT_path_by_regexp_no_match(testCase)
      bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
        testCase, {'20201.cdf', '2020.CDF'});

      testCase.verifyError(...
        @() bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        testCase.dir, '20[0-9][0-9]\.cdf', testCase.L), ...
        'BICAS:CannotFindRegexMatchingRCT')
    end



    function test_get_RCT_path_by_regexp_1_match(testCase)
      bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
        testCase, {'2020.cdf', 'asdsf'});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        testCase.dir, '20[0-9][0-9]\.cdf', testCase.L);

      testCase.verifyEqual(...
        path, ...
        fullfile(testCase.dir, '2020.cdf'))
    end



    function test_get_RCT_path_by_regexp_2_match(testCase)
      bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
        testCase, {'2020.cdf', '2021.cdf'});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        testCase.dir, '20[0-9][0-9]\.cdf', testCase.L);

      testCase.verifyEqual(...
        path, ...
        fullfile(testCase.dir, '2021.cdf'))
    end



    function test_get_RCT_path_by_regexp_realistic(testCase)
      FN_1 = 'SOLO_CAL_RPW-BIAS_V202111191204.cdf';
      FN_2 = 'SOLO_CAL_RPW-BIAS_V202011191204.cdf';

      bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
        testCase, {FN_1, FN_2});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        testCase.dir, 'SOLO_CAL_RPW-BIAS_V20[0-9]{10,10}.cdf', testCase.L);

      testCase.verifyEqual(...
        path, ...
        fullfile(testCase.dir, FN_1))
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function setup_files(testCase, filenamesCa)
      for fileCa = filenamesCa(:)'
        irf.fs.write_empty_file({testCase.dir, fileCa{1}});
      end
    end



  end    % methods(Static, Access=private)



end
