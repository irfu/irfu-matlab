%
% matlab.unittest automatic test code for bicas.proc.L1L2.SdChannelsData.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SdChannelData___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor_size_bFp(testCase)
      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(0, 0), false(0, 1));
      testCase.assertEqual(size(Sdcd), [0, 1])
      testCase.assertEqual(Sdcd.bFp,   false(0, 1))

      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(0, 3), false(0, 1));
      testCase.assertEqual(size(Sdcd), [0, 1])
      testCase.assertEqual(Sdcd.bFp,   false(0, 1))

      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(1, 3), false);
      testCase.assertEqual(size(Sdcd), [1, 1])
      testCase.assertEqual(Sdcd.bFp,   false(1, 1))



      SAMPLES_AR = [1,2; 3,4; 5,NaN];
      VSQB_AR    = logical([0;1;0]);

      Sdcd = bicas.proc.L1L2.SdChannelData(SAMPLES_AR, VSQB_AR);
      testCase.assertEqual(Sdcd.samplesAr, SAMPLES_AR)
      testCase.assertEqual(Sdcd.vsqbAr,    VSQB_AR)
      testCase.assertEqual(size(Sdcd),     [3, 1])
      testCase.assertEqual(Sdcd.bFp,       logical([0; 0; 1]))
    end



    function test_subsref(testCase)
      SDCD_123 = testCase.get_SDCD([1,2,3]);
      SDCD_13  = testCase.get_SDCD([1,3]);
      SDCD_2   = testCase.get_SDCD([2]);
      SDCD_    = testCase.get_SDCD([]);

      testCase.assertEqual(SDCD_123, SDCD_123(1:3))
      testCase.assertEqual(SDCD_123, SDCD_123(logical([1,1,1])))

      testCase.assertEqual(SDCD_,  SDCD_([]))
      testCase.assertEqual(SDCD_,  SDCD_(logical([])))

      testCase.assertEqual(SDCD_13, SDCD_123(logical([1,0,1])))
      testCase.assertEqual(SDCD_13, SDCD_123([1,3]))

      testCase.assertEqual(SDCD_2,  SDCD_123(logical([0,1,0])))
      testCase.assertEqual(SDCD_2,  SDCD_123([2]))

      testCase.assertEqual(SDCD_,  SDCD_123(logical([0,0,0])))
      testCase.assertEqual(SDCD_,  SDCD_123([]))

      % % BELOW SHOULD BE EQUIVALENT TO PARTS OF CODE ABOVE. TOO CONVOLUTED?
      % function test(ibExp, ibSdcdInit, ibSubsref)
      %   ExpSdcd = testCase.get_SDCD(ibExp);
      %   Sdcd    = testCase.get_SDCD(ibSdcdInit);
      %   ActSdcd = Sdcd(ibSubsref);
      %   testCase.assertEqual(ActSdcd, ExpSdcd)
      % end
      %
      % test(1:3, 1:3, 1:3)
      % test(1:3, 1:3, logical([1,1,1]))
      %
      % test([], [], [])
      % test([], [], logical([]))
    end



    function test_subasgn(testCase)
      function test()
        SDCD_R          = testCase.get_SDCD(ib_r);
        SDCD_S          = testCase.get_SDCD(ib_s);
        SDCD_EXP        = testCase.get_SDCD(ib_exp);
        SDCD_R(ib_asgn) = SDCD_S;
        testCase.assertEqual(SDCD_R, SDCD_EXP)
      end

      % Size 0-->1
      ib_r    = [];
      ib_s    = [3];
      ib_asgn = [true];
      ib_exp  = [3];
      test()
      ib_r    = [];
      ib_s    = [3];
      ib_asgn = [1];
      ib_exp  = [3];
      test()

      % Overwrite 1 of 3.
      ib_r    = [];
      ib_s    = [3];
      ib_asgn = [true];
      ib_exp  = [3];
      test()
      % Overwrite 1 of 3, add 1.
      ib_r    = [1,2,3];
      ib_s    = [4,5];
      ib_asgn = [2,4];
      ib_exp  = [1,4,3,5];
      test()
      ib_r    = [1,2,3];
      ib_s    = [4,5];
      ib_asgn = logical([0,1,0,1]);
      ib_exp  = [1,4,3,5];
      test()
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Get test SDCD with test data by indexing input to SDCD.
    function Sdcd = get_SDCD(testCase, ib)
      SAMPLES_AR = [1,2; 3,4; 5,6; 7,8; 9,10];
      VSQB_AR    = logical([0;1;0;1;0]);

      Sdcd = bicas.proc.L1L2.SdChannelData(SAMPLES_AR(ib, :), VSQB_AR(ib, :));
    end



  end    % methods(Access=private)



end
