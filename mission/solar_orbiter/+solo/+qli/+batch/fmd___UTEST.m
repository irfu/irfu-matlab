%
% matlab.unittest automatic test code for solo.qli.batch.fmd.
%
% NOTE: Does not test (very much) supported functionality that is not expected
% to be really needed.
% * Non-24 h QLI filenames.
% * Non-24 h datasets (by filename).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef fmd___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_days_from_QLI_FMD_interval(testCase)

      function test(FsrDict, qliDir, startInclFmdStr, stopExclFmdStr, ExpQliDfmdd)
        Fsr = solo.qli.batch.FileSystemReaderTest(FsrDict);

        ActQliDfmdd = solo.qli.batch.fmd.get_days_from_QLI_FMD_interval(...
          qliDir, datetime(startInclFmdStr), datetime(stopExclFmdStr), Fsr);

        testCase.assertEqual(ActQliDfmdd, ExpQliDfmdd)
      end

      %=========================================================================

      if 1
        % Zero QLIs.
        FsrDict = dictionary();
        FsrDict({{'/qli'}}) = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
        ExpDfmdd = dictionary();
        test( ...
          FsrDict, ...
          '/qli', ...
          '2020-01-01', '2020-02-01', ...
          ExpDfmdd ...
          )
      end



      % Non-zero QLIs
      FsrDict = dictionary();
      QLI_CA = {
        '20240101T00_20240102T00.png'; ...
        '20240102T00_20240103T00.png';
        '20240103T00_20240104T00.png';
        '20240104T00_20240105T00.png';
        };
      FMD_DT_ARRAY = datetime({
        '2025-01-01';
        '2025-01-02';
        '2025-01-03';
        '2025-01-04';
        });
      FsrDict({{'/qli'}}) = {{QLI_CA, FMD_DT_ARRAY}};

      ExpDfmdd = dictionary();
      ExpDfmdd(solo.qli.utils.umdt({ ...
        '2024-01-02';
        '2024-01-03'
        })) = datetime({ ...
        '2025-01-02';
        '2025-01-03'
      });

      test( ...
        FsrDict, ...
        '/qli', ...
        '2025-01-02', '2025-01-04', ...
        ExpDfmdd ...
        )

    end



    % Test using actual file system (FS). ==> More complicated by nature.
    function test_get_days_from_IDMRQ_and_FS_FS(testCase)
      QliFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      fmdQliDir  = QliFixture.Folder;
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir1 = F.Folder;
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir2 = F.Folder;

      % NOTE: Do NOT use test class.
      Fsr = solo.qli.batch.FileSystemReaderImplementation();



      %================================
      % Create quicklooks and datasets
      %================================
      % NOTE: Quicklook (created first) is older than datasets.
      [~] = irf.fs.create_empty_file({fmdQliDir, '20240101T00_20240102T00.png'});

      % Delay does not appear to be needed. Added for safety.
      pause(0.1)

      [~] = irf.fs.create_empty_file({dir1, 'solo_L2_swa-pas-eflux_20240101_V02.cdf'});
      [~] = irf.fs.create_empty_file({dir2, 'solo_L2_mag-rtn-normal_20240102_V02.cdf'});

      % Should be needed since FMDs appear to have a time-resolution of 1 s (?)
      % (Linux).
      % NOTE: This delay is a potential source of mistaken test failures if the
      % local file modification time or clock does not support a time resolution
      % smaller than this delay.
      pause(1.1)

      % NOTE: Quicklook (created last) is more recent than datasets.
      [~] = irf.fs.create_empty_file({fmdQliDir, '20240102T00_20240103T00.png'});



      %==============
      % Execute test
      %==============
      datasetDirsCa = {dir1; dir2};

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_IDMRQ(...
        datasetDirsCa, fmdQliDir, Fsr, ...
        {'9999', '2000-01-01', '2099-01-01'});

      testCase.assertEqual(ActDaysDtArray, ...
        solo.qli.utils.umdt({'2024-01-01'}))
    end



    % Test using solo.qli.batch.FileSystemReaderTest. ==> Less complicated by
    % nature.
    function test_get_days_from_IDMRQ_and_FS_FSRTest(testCase)
      function test(datasetDirsCa, qliDir, dsiCa, FsrDict, expDaysStrCa)
        assert(ischar(qliDir))

        Fsr            = solo.qli.batch.FileSystemReaderTest(FsrDict);
        ExpDaysDtArray = solo.qli.utils.umdt(expDaysStrCa);

        ActDaysDtArray = solo.qli.batch.fmd.get_days_from_IDMRQ_and_FS(...
          datasetDirsCa, qliDir, dsiCa, Fsr);

        testCase.assertEqual(ActDaysDtArray, ExpDaysDtArray)
      end

      %=========================================================================

      % Empty
      FsrDict = dictionary();
      FsrDict({cell(0,1)}) = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      FsrDict({{'/qli'}})  = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      test(...
        cell(0, 1), ...
        '/qli', ...
        cell(0, 1), ...
        FsrDict, ...
        cell(0, 1)...
        )

      % One dataset, no QLI
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{'/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'}, datetime('2025-01-01')}};
      FsrDict({{'/qli' }}) = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      test( ...
        {'/data'}, ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX'}, ...
        FsrDict, ...
        cell(0, 1) ...
        )

      % No dataset, one QLI
      FsrDict = dictionary();
      FsrDict({cell(0, 1)}) = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      FsrDict({{'/qli'}})   = {{{'/data/20240101T00_20240102T00.png'}, datetime('2025-01-01')}};
      test( ...
        cell(0, 1), ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX'}, ...
        FsrDict, ...
        cell(0, 1) ...
        )

      % One dataset younger than QLI.
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{'/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'}, datetime('2025-01-02')}};
      FsrDict({{'/qli' }}) = {{{'/qli/20240101T00_20240102T00.png'            }, datetime('2025-01-01')}};
      test( ...
        {'/data'}, ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX'}, ...
        FsrDict, ...
        {'2024-01-01'} ...
        )

      % One dataset older than QLI.
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{'/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'}, datetime('2025-01-01')}};
      FsrDict({{'/qli' }}) = {{{'/qli/20240101T00_20240102T00.png'            }, datetime('2025-01-02')}};
      test( ...
        {'/data'}, ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX'}, ...
        FsrDict, ...
        cell(0, 1) ...
        )

      % One day, two datasets, two DSIs, one younger and one older than QLI.
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{...
        '/data/solo_L2_swa-pas-eflux_20240101_V02.cdf';
        '/data/solo_L2_mag-rtn-normal_20240101_V02.cdf' ...
        }, ...
        datetime({ ...
        '2025-01-01'; ...
        '2025-01-03' ...
        })}};
      FsrDict({{'/qli' }}) = {{{'/qli/20240101T00_20240102T00.png'}, datetime('2025-01-02')}};
      test( ...
        {'/data'}, ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX', 'SOLO_L2_MAG-RTN-NORMAL'}, ...
        FsrDict, ...
        {'2024-01-01'} ...
        )

      % One day, two datasets, two DSIs, both older than QLI.
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{...
        '/data/solo_L2_swa-pas-eflux_20240101_V02.cdf';
        '/data/solo_L2_mag-rtn-normal_20240101_V02.cdf' ...
        },
        datetime({ ...
        '2025-01-03'; ...
        '2025-01-01' ...
        })}};
      FsrDict({{'/qli' }}) = {{{'/qli/20240101T00_20240102T00.png'}, datetime('2025-01-04')}};
      test( ...
        {'/data'}, ...
        '/qli', ...
        {'SOLO_L2_SWA-PAS-EFLUX', 'SOLO_L2_MAG-RTN-NORMAL'}, ...
        FsrDict, ...
        cell(0, 1) ...
        )
    end



    % NOTE: Not solo.qli.batch.fmd.testing get_dataset_DFMDD_for_one_DSI()
    % which is wrapped by
    % solo.qli.batch.fmd.get_dataset_DFMDD_for_all_DSIs().
    function test_get_dataset_DFMDD_for_all_DSIs(testCase)

      function test(datasetPathsCa, fmdStrCa, dsiCa, expDaysStrCa, expFmdStrCa)
        [DsmdArray, ~] = solo.adm.paths_to_DSMD_array(datasetPathsCa);
        FmdDtArray = datetime(fmdStrCa);
        ExpDayFmdDict = dictionary(...
          solo.qli.utils.umdt(expDaysStrCa), ...
          datetime(           expFmdStrCa));

        ActDayFmdDict = solo.qli.batch.fmd.get_dataset_DFMDD_for_all_DSIs(...
          DsmdArray, FmdDtArray, dsiCa);

        testCase.assertEqual(ActDayFmdDict, ExpDayFmdDict)
      end

      %=========================================================================

      % Empty
      test(...
        cell(0, 1), ...
        cell(0, 1), ...
        cell(0, 1), ...
        cell(0, 1), ...
        cell(0, 1))

      % File does not match specified DSI.
      % DSI does not match file.
      test(...
        {'/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'}, ...
        {'2024-01-01'}, ...
        {'SOLO_L2_MAG-RTN-NORMAL'}, ...
        cell(0, 1), ...
        cell(0, 1))

      % Two datasets: same day, different DSIs and FMDs.
      test(...
        { ...
        '/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'; ...
        '/data/solo_L2_mag-rtn-normal_20240101_V02.cdf'; ...
        }, ...
        { ...
        '2025-01-02';
        '2025-01-01';
        }, ...
        {'SOLO_L2_SWA-PAS-EFLUX', 'SOLO_L2_MAG-RTN-NORMAL'}, ...
        {'2024-01-01'}, ...
        {'2025-01-02'} ...
        )
    end



    function test_construct_QLI_DFMDD(testCase)

      function test(qliPathsCa, QliFmdDtArray, ExpDayFmdDict)
        ActDayFmdDict = solo.qli.batch.fmd.construct_QLI_DFMDD(...
          qliPathsCa, QliFmdDtArray);

        testCase.assertEqual(ActDayFmdDict, ExpDayFmdDict)
      end

      %=========================================================================

      test( ...
        cell(0, 1), ...
        NaT(0, 1), ...
        dictionary(datetime.empty, datetime.empty) ...
        )

      % One non-matching filename
      test( ...
        {'illegal_quicklook_filename'}, ...
        datetime('27-Mar-2024 14:27:39'), ...
        dictionary(datetime.empty, datetime.empty) ...
        )

      % One matching filename: Day-long, midnight-midnight.
      FMD_DT = datetime('27-Mar-2024 14:27:39');
      test( ...
        {'20220101T00_20220102T00.png'}, ...
        FMD_DT, ...
        dictionary(solo.qli.utils.umdt('2022-01-01'), FMD_DT) ...
        )

      % One matching filename: Non-midnight to non-midnight.
      FMD_DT = datetime('27-Mar-2024 14:27:39');
      test( ...
        {'20220101T22_20220102T02.png'}, ...
        FMD_DT, ...
        dictionary(...
        solo.qli.utils.umdt({'2022-01-01'; '2022-01-02'}), ...
        [FMD_DT; FMD_DT]) ...
        )

      % One matching filename: Two-day-long, midnight-midnight.
      FMD_DT = datetime('27-Mar-2024 14:27:39');
      test( ...
        {'20220101T00_20220103T00.png'}, ...
        FMD_DT, ...
        dictionary(...
        solo.qli.utils.umdt({'2022-01-01'; '2022-01-02'}), ...
        FMD_DT))

      % "Complex test"
      FMD_DT_ARRAY = datetime({ ...
        '27-Mar-2024 14:00:00';
        '27-Mar-2024 15:00:00';
        '27-Mar-2024 16:00:00';
        });
      test( ...
        { ...
        '20220101T00_20220103T00.png';   % 2 days
        '20220102T00_20220103T00.png';   % 1 day that overlaps.
        '20220103T00_20220104T00.png';   % 1 day that overlaps with nothing.
        }, ...
        FMD_DT_ARRAY, ...
        dictionary(...
        solo.qli.utils.umdt({ ...
        '2022-01-01'; ...
        '2022-01-02'; ...
        '2022-01-03'}), ...
        FMD_DT_ARRAY) ...
        )
    end



  end    % methods(Test)



end
