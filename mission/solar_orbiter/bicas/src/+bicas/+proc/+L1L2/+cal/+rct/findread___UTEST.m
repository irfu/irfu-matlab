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



    function setup(T)
      Fixture = T.applyFixture(...
        matlab.unittest.fixtures.TemporaryFolderFixture);
      % NOTE: The same fixture should always return the same directory.
      T.dir = Fixture.Folder;
      T.L   = bicas.Logger('NO_STDOUT', false);
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_BRVF_RCT_path(T)
      % Create BRVF.
      ExpDtBegin = datetime('2020-02-10T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDtEnd   = datetime('2100-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      expBiasRctFilename = bicas.tools.rct.create_RCT_filename(ExpDtBegin, ExpDtEnd, 1);
      bicas.tools.rct.create_BRVF(T.dir, expBiasRctFilename, ExpDtBegin, ExpDtEnd);

      expBiasRctPath = fullfile(T.dir, expBiasRctFilename);
      expBrvfPath    = fullfile(T.dir, bicas.const.BRVF_FILENAME);

      %========
      % TEST 1
      %========
      [actBiasRctPath, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.get_BRVF_RCT_path(T.dir, ExpDtBegin, ExpDtEnd);

      T.assertEqual(actBiasRctPath, expBiasRctPath)
      T.assertEqual(actBrvfPath,    expBrvfPath)



      %========
      % TEST 2
      %========
      Duration = duration('00:00:10');
      % NOTE: Only checks for non-error time interval case.
      [actBiasRctPath, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.get_BRVF_RCT_path(...
        T.dir, ...
        ExpDtBegin + Duration, ...
        ExpDtEnd   - Duration);
    end



    function test_read_BRVF___compliant_RCT_filename(T)
      % Create BRVF.
      ExpDtBegin = datetime('2020-02-10T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDtEnd   = datetime('2100-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      expBiasRctFilename = bicas.tools.rct.create_RCT_filename(ExpDtBegin, ExpDtEnd, 1);
      bicas.tools.rct.create_BRVF(T.dir, expBiasRctFilename, ExpDtBegin, ExpDtEnd);

      [actRctFilename, ActDtValidityBegin, ActDtValidityEnd, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.read_BRVF(T.dir);

      T.assertEqual(actRctFilename,     expBiasRctFilename)
      T.assertEqual(ActDtValidityBegin, ExpDtBegin)
      T.assertEqual(ActDtValidityEnd,   ExpDtEnd)
      T.assertEqual(actBrvfPath,        fullfile(T.dir, bicas.const.BRVF_FILENAME))
    end



    function test_read_BRVF___noncompliant_RCT_filename(T)
      % Create BRVF.
      ExpDtBegin = datetime('2023-04-05T06:07:08Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDtEnd   = datetime('2100-02-03T04:05:06Z', 'TimeZone', 'UTCLeapSeconds');
      expBiasRctFilename = 'SOLO_CAL_RPW-BIAS_V202011191204.cdf';   % Non-compliant RCT filename.
      bicas.tools.rct.create_BRVF(T.dir, expBiasRctFilename, ExpDtBegin, ExpDtEnd);

      [actRctFilename, ActDtValidityBegin, ActDtValidityEnd, actBrvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.read_BRVF(T.dir);

      T.assertEqual(actRctFilename,     expBiasRctFilename)
      T.assertEqual(ActDtValidityBegin, ExpDtBegin)
      T.assertEqual(ActDtValidityEnd,   ExpDtEnd)
      T.assertEqual(actBrvfPath,        fullfile(T.dir, bicas.const.BRVF_FILENAME))
    end



    function test_get_RCT_path_by_regexp_empty(T)
      T.setup_files({});

      T.verifyError(...
        @() bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        T.dir, '20[0-9][0-9]\.cdf', T.L), ...
        'BICAS:CannotFindRegexMatchingRCT')

    end



    function test_get_RCT_path_by_regexp_no_match(T)
      T.setup_files(...
        {'20201.cdf', '2020.CDF'});

      T.verifyError(...
        @() bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        T.dir, '20[0-9][0-9]\.cdf', T.L), ...
        'BICAS:CannotFindRegexMatchingRCT')
    end



    function test_get_RCT_path_by_regexp_1_match(T)
      T.setup_files(...
        {'2020.cdf', 'asdsf'});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        T.dir, '20[0-9][0-9]\.cdf', T.L);

      T.verifyEqual(...
        path, ...
        fullfile(T.dir, '2020.cdf'))
    end



    function test_get_RCT_path_by_regexp_2_match(T)
      T.setup_files(...
        {'2020.cdf', '2021.cdf'});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        T.dir, '20[0-9][0-9]\.cdf', T.L);

      T.verifyEqual(...
        path, ...
        fullfile(T.dir, '2021.cdf'))
    end



    function test_get_RCT_path_by_regexp_realistic(T)
      FN_1 = 'SOLO_CAL_RPW-BIAS_V202111191204.cdf';
      FN_2 = 'SOLO_CAL_RPW-BIAS_V202011191204.cdf';

      T.setup_files(...
        {FN_1, FN_2});

      path = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        T.dir, 'SOLO_CAL_RPW-BIAS_V20[0-9]{10,10}.cdf', T.L);

      T.verifyEqual(...
        path, ...
        fullfile(T.dir, FN_1))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function setup_files(T, filenamesCa)
      for fileCa = filenamesCa(:)'
        irf.fs.write_empty_file({T.dir, fileCa{1}});
      end
    end



  end    % methods(Access=private)



end
