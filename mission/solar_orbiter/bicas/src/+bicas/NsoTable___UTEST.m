%
% matlab.unittest automatic test code for bicas.NsoTable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef NsoTable___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Extend with more tests.
  %   PROPOSAL: Empty (legal) NSO file?
  %       NOTE: Already test-loading default NSO XML file.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % % Test constructor.
    % % NOTE: test_get_NSO_events_timestamps() indirectly tests the constructor.
    % function test_NsoTable(testCase)
    % end



    % Test method bicas.NsoTable.get_NSO_events_timestamps(). Also indirectly
    % uses/tests bicas.NsoTable() constructor by its nature.
    function test_get_NSO_events_timestamps(testCase)
      % PROBLEM: How handle that return value may change the order of
      %          events depending on implementation?

      function test(...
          evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr, tt2000Array, ...
          expBEvtArraysCa, expEvtQrcidAr, expIGlobalEventsArray)

        % Normalize input
        evtStartTt2000Array = int64(evtStartTt2000Array(:));
        evtStopTt2000Array  = int64(evtStopTt2000Array(:));
        evtQrcidAr          = evtQrcidAr(:);
        tt2000Array         = int64(tt2000Array(:));
        % Normalize (expected) output
        expBEvtArraysCa       = expBEvtArraysCa(:);
        expEvtQrcidAr         = expEvtQrcidAr(:);
        expIGlobalEventsArray = expIGlobalEventsArray(:);

        NsoTable = bicas.NsoTable(...
          evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr);

        [actBEvtArraysCa, actEvtQrcidAr, actIGlobalEventsArray] = ...
          NsoTable.get_NSO_events_timestamps(tt2000Array);
        testCase.verifyEqual(actBEvtArraysCa,       expBEvtArraysCa)
        testCase.verifyEqual(actEvtQrcidAr,         expEvtQrcidAr)
        testCase.verifyEqual(actIGlobalEventsArray, expIGlobalEventsArray)
      end

      %===================================================================

      QRCID_1 = "QRCID_1_FOR_TESTING";
      QRCID_2 = "QRCID_2_FOR_TESTING_MORE";

      ESA = strings(0, 1);   % ESA=Empty String  Array
      ENA = zeros(  0, 1);   % ENA=Empty Numeric Array
      ECA = cell(   0, 1);   % ECA=Empty Cell    Array

      % Test empty output
      % -----------------
      % Test every combination of
      % (1) empty & (2) non-empty
      % for (a) NSO table & (b) submitted timestamps
      % without any overlap (empty output).
      for tt2000ArrayCa = {ENA, [100:200]'}
        tt2000Ar = tt2000ArrayCa{1};

        % Empty NSO table.
        test(...
          ENA, ENA, ESA, ...
          tt2000Ar, ...
          ECA, ESA, ENA)

        % Non-empty NSO table.
        test(...
          [10], [20], QRCID_1, ...
          tt2000Ar, ...
          ECA, ESA, ENA)
      end

      % NSO events do not overlap beginning & end.
      test(...
        [1, 5], [2, 7], [QRCID_1, QRCID_2], ...
        [0:9], ...
        {...
        %        0  1  2  3  4  5  6  7  8  9
        logical([0, 1, 1, 0, 0, 0, 0, 0, 0, 0]'), ...
        logical([0, 0, 0, 0, 0, 1, 1, 1, 0, 0]') ...
        }, ...
        [QRCID_1, QRCID_2], [1, 2])

      % NSO events overlap beginning & end.
      test(...
        [-1, 5], [2, 12], [QRCID_1, QRCID_2], ...
        [0:9], ...
        {...
        %        0  1  2  3  4  5  6  7  8  9
        logical([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]'), ...
        logical([0, 0, 0, 0, 0, 1, 1, 1, 1, 1]') ...
        }, ...
        [QRCID_1, QRCID_2], [1, 2])

      % Overlapping NSOs.
      test(...
        [ 1, 3], [5, 8], [QRCID_1, QRCID_2], ...
        [0:9], ...
        {...
        %        0  1  2  3  4  5  6  7  8  9
        logical([0, 1, 1, 1, 1, 1, 0, 0, 0, 0]'), ...
        logical([0, 0, 0, 1, 1, 1, 1, 1, 1, 0]') ...
        }, ...
        [QRCID_1, QRCID_2], [1, 2])

    end



    function test_read_file_validated(testCase)
      % NOTE: Only read BICAS's own default file (in BICAS's git repo).

      bicasRootPath = bicas.utils.get_BICAS_root_dir();

      Bso = bicas.create_default_BSO();
      Bso.make_read_only()
      nsoTableRelativePath = Bso.get_fv('PROCESSING.NSO_TABLE.FILE.RELATIVE_PATH');

      nsoFilePath = fullfile(bicasRootPath, nsoTableRelativePath);

      % CALL TESTED CODE
      NsoTable = bicas.NsoTable.read_file_validated(nsoFilePath, bicas.const.ALL_QRCID_AR);

      testCase.verifyTrue(isa(NsoTable, 'bicas.NsoTable'))
      testCase.assertTrue(isstring(NsoTable.evtQrcidAr))

      nEvents = irf.assert.sizes( ...
        NsoTable.evtStartTt2000Array, [-1, 1], ...
        NsoTable.evtStopTt2000Array,  [-1, 1], ...
        NsoTable.evtQrcidAr,          [-1, 1]);
      testCase.verifyTrue(nEvents > 300)
    end



  end    % methods(Test)



end
