%
% matlab.unittest automatic test code for bicas.proc.L1L2.demuxer.
%
% Could be improved but unsure how much is meaningful. Seems to complicated.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, using older test code.
%
classdef demuxer___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_reconstruct_missing_data(testCase)

      % Test data with relationship a1 + a2 == a3 for non-NaN.
      function test_sum(ACa, expACa)
        [actA1,actA2,actA3] = bicas.proc.L1L2.demuxer.reconstruct_missing_data( ...
          ACa{1},        ACa{2},        ACa{3}, ...
          isnan(ACa{1}), isnan(ACa{2}), isnan(ACa{3}), ...
          @(a2, a3) (a3-a2), ...
          @(a1, a3) (a3-a1), ...
          @(a1, a2) (a1+a2));

        testCase.assertEqual(actA1, expACa{1})
        testCase.assertEqual(actA2, expACa{2})
        testCase.assertEqual(actA3, expACa{3})
      end

      %==============
      % Empty arrays
      %==============
      for sizeCa = {[0,0], [1,0], [0,1], [0,0,1]}
        A = zeros(sizeCa{1});
        test_sum({A, A, A}, {A, A, A})
      end

      %====================
      % Scalars, all cases
      %====================
      % Inconsistent existing relationship.
      test_sum({1, 3, 9}, {1, 3, 9})

      % Reconstruct
      test_sum({nan, 3, 5}, {2, 3, 5})
      test_sum({2, nan, 5}, {2, 3, 5})
      test_sum({2, 3, nan}, {2, 3, 5})

      % Can not reconstruct values.
      test_sum({2, nan, nan}, {2, nan, nan})
      test_sum({nan, 3, nan}, {nan, 3, nan})
      test_sum({nan, nan, 5}, {nan, nan, 5})
      %
      test_sum({nan, nan, nan}, {nan, nan, nan})

      %==========================================
      % Non-scalar array, multiple cases at once
      %==========================================
      test_sum({ ...
        [1, nan,   4;   7, nan, nan], ...
        [3,   3, nan;   8,   9, nan], ...
        [9,   5,  10; nan, nan, nan] ...
        }, { ...
        [1,   2,   4;   7, nan, nan], ...
        [3,   3,   6;   8,   9, nan], ...
        [9,   5,  10;  15, nan, nan] ...
        })
    end



    function test_reconstruct_ASR_samples_NEW(testCase)
      N = NaN;

      SAMPLES_AR_DATA = [...
        1,2,4,   -1,-3,-2,   7,10,3;
        1,N,N,   -1, N,-2,   7,10,N;
        1,N,N,   -1, N,-2,   7, N,N;
        ];
      VSIB_AR_DATA = logical([ ...
        1,0,0,    0, 0, 0,   1, 0,0;
        1,0,0,    0, 0, 0,   1, 0,0;
        0,0,0,    0, 0, 1,   1, 0,0;
        ]);
      Zvm = testCase.create_ASR_SCDH_ZVM(SAMPLES_AR_DATA, VSIB_AR_DATA);



      EXP_SAMPLES_AR_DATA = [...
        1,2,4,   -1,-3,-2,   7,10,3;
        1,2,4,   -1,-3,-2,   7,10,3;
        1,2,4,   -1,-3,-2,   7, N,N;
        ];
      EXP_VSIB_AR_DATA = logical([ ...
        1,0,0,    0, 0, 0,   1, 0,0;
        1,1,1,    0, 0, 0,   1, 0,1;
        0,0,1,    0, 1, 1,   1, 0,0;
        ]);
      ExpZvm = testCase.create_ASR_SCDH_ZVM(EXP_SAMPLES_AR_DATA, EXP_VSIB_AR_DATA);



      bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_NEW(Zvm);
      ActZvm = Zvm;

      testCase.assertEqual(ActZvm.nEntries, 9)

      % IMPLEMENTATION NOTE: Not only comparing entire ZVM objects since it
      % helps to compare object components separately when debugging tests.
      for sdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'
        ActSchd = ActZvm.get(sdid);
        ExpSchd = ExpZvm.get(sdid);

        % Print/log component values if not equal (for debugging).
        if ~isequaln(ActSchd.samplesAr, ExpSchd.samplesAr)
          sdid
          ActSchd.samplesAr
          ExpSchd.samplesAr
        end
        if ~isequaln(ActSchd.vsibAr, ExpSchd.vsibAr)
          sdid
          ActSchd.vsibAr
          ExpSchd.vsibAr
        end

        % Check everything (partially overlapping with above).
        testCase.assertEqual(ActSchd, ExpSchd)
      end
      testCase.assertEqual(ActZvm, ExpZvm)
    end



    function test_get_ASR_ZVM_nWholeRowIsNan(testCase)
      SDID_ASR_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR;

      % 3x2
      SAMPLES_AR_0 = [NaN,2; 3,NaN; NaN,NaN];
      VSIB_AR      = logical([0; 1; 0]);

      Zvm = bicas.utils.ZvMap(3);

      for i = 1:numel(SDID_ASR_AR)
        sdid = SDID_ASR_AR(i);

        Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR_0+i, VSIB_AR);
        Zvm.add(sdid, Schd);

        actNWholeRowIsNan = bicas.proc.L1L2.demuxer.get_ASR_ZVM_nWholeRowIsNan(Zvm);
        testCase.assertEqual(actNWholeRowIsNan, 1*i)
      end

      % Double checks. Not needed for the test.
      testCase.assertEqual(Zvm.nRecords, 3)
      ExpSchd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR_0+3, VSIB_AR);
      ActSchd = Zvm.get(SDID_ASR_AR(3));
      testCase.assertEqual(ActSchd, ExpSchd)

      actNWholeRowIsNan = bicas.proc.L1L2.demuxer.get_ASR_ZVM_nWholeRowIsNan(Zvm);
      testCase.assertEqual(actNWholeRowIsNan, numel(SDID_ASR_AR)*1)
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Fast-and-easy function for creating one ZVM SDID-->SCHD from variables on
    % a format suitable for hardcoding (CWF only).
    %
    % ARGUMENTS
    % =========
    % samplesArData
    %       (iRec, iSdid). Only column arrays will be assigned to each separate
    %       SCHD.
    % vsibArData
    %       (iRec, iSdid)
    %
    function Zvm = create_ASR_SCDH_ZVM(samplesArData, vsibArData)
      assert(size(samplesArData, 2) == 9)
      assert(size(vsibArData,    2) == 9)

      SDID_AR  = bicas.proc.L1L2.const.C.SDID_ASR_AR;
      nRecords = size(vsibArData, 1);

      Zvm = bicas.utils.ZvMap(nRecords);
      for iSdid = 1:numel(SDID_AR)

        Schd = bicas.proc.L1L2.SingleChannelData(...
          samplesArData(:, iSdid), ...
          vsibArData(   :, iSdid));
        Zvm.add(SDID_AR(iSdid), Schd);
      end
    end



  end



end
