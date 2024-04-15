%
% matlab.unittest automatic test code for
% solo.qli.batch.interface.get_days_from_DMRQ().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface_get_days_from_DMRQ___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(testCase)
      FsrDict = dictionary();
      FsrDict({{'/qli'}})  = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      FsrDict({cell(0,1)}) = {{cell(0, 1), solo.qli.const.EMPTY_DT_ARRAY}};
      Fsr = solo.qli.batch.FileSystemReaderTest(FsrDict);

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_DMRQ(...
        cell(0, 1), '/qli', Fsr, ...
        {'999', '2000-01-01', '2099-01-01'});

      testCase.assertEqual(ActDaysDtArray, solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_nonempty(testCase)
      % One dataset is younger than QLI.
      FsrDict = dictionary();
      FsrDict({{'/data'}}) = {{{'/data/solo_L2_swa-pas-eflux_20240101_V02.cdf'}, datetime('2025-01-02')}};
      FsrDict({{'/qli' }}) = {{{'/qli/20240101T00_20240102T00.png'            }, datetime('2025-01-01')}};
      Fsr = solo.qli.batch.FileSystemReaderTest(FsrDict);

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_DMRQ(...
        {'/data'}, '/qli', Fsr, ...
        {'999', '2000-01-01', '2099-01-01'});

      testCase.assertEqual(ActDaysDtArray, solo.qli.utils.umddt({'2024-01-01'}))
    end



  end    % methods(Test)



end
