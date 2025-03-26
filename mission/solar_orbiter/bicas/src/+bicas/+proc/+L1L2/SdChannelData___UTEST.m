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



    function test_constructor_size_bWholeRowIsNan(testCase)
      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(0, 0), false(0, 1));
      testCase.assertEqual(size(Sdcd),  [0, 1])
      testCase.assertEqual(Sdcd.bWholeRowIsNan, false(0, 1))

      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(0, 3), false(0, 1));
      testCase.assertEqual(size(Sdcd),  [0, 1])
      testCase.assertEqual(Sdcd.bWholeRowIsNan, false(0, 1))

      Sdcd = bicas.proc.L1L2.SdChannelData(zeros(1, 3), false(1, 1));
      testCase.assertEqual(size(Sdcd),  [1, 1])
      testCase.assertEqual(Sdcd.bWholeRowIsNan, false(1, 1))



      SAMPLES_AR = [1,2; 3,4; 5,NaN; NaN,8; NaN,NaN];
      VSQB_AR    = logical([0; 1; 0; 1; 0]);

      Sdcd = bicas.proc.L1L2.SdChannelData(SAMPLES_AR, VSQB_AR);
      testCase.assertEqual(Sdcd.samplesAr,      SAMPLES_AR)
      testCase.assertEqual(Sdcd.vsqbAr,         VSQB_AR)
      testCase.assertEqual(size(Sdcd),          [5, 1])
      testCase.assertEqual(Sdcd.bWholeRowIsNan, logical([0; 0; 0; 0; 1]))
    end



    function test_subsref(testCase)
      SDCD_123 = testCase.get_SDCD([1,2,3]);
      SDCD_13  = testCase.get_SDCD([1,3]);
      SDCD_31  = testCase.get_SDCD([3,1]);
      SDCD_2   = testCase.get_SDCD([2]);
      SDCD_    = testCase.get_SDCD([]);

      testCase.assertEqual(SDCD_123, SDCD_123(1:3))
      testCase.assertEqual(SDCD_123, SDCD_123(logical([1,1,1])))

      testCase.assertEqual(SDCD_,   SDCD_([]))
      testCase.assertEqual(SDCD_,   SDCD_(logical([])))

      testCase.assertEqual(SDCD_13, SDCD_123(logical([1,0,1])))
      testCase.assertEqual(SDCD_13, SDCD_123([1,3]))

      testCase.assertEqual(SDCD_31, SDCD_123([3,1]))

      testCase.assertEqual(SDCD_2,  SDCD_123(logical([0,1,0])))
      testCase.assertEqual(SDCD_2,  SDCD_123([2]))

      testCase.assertEqual(SDCD_,   SDCD_123(logical([0,0,0])))
      testCase.assertEqual(SDCD_,   SDCD_123([]))
    end



    function test_subasgn(testCase)
      function test()
        % R = Receiver?!
        % S = Sender?!
        SDCD_R          = testCase.get_SDCD(ib_r);
        SDCD_S          = testCase.get_SDCD(ib_s);
        SDCD_EXP        = testCase.get_SDCD(ib_exp);
        SDCD_R(ib_r_asgn) = SDCD_S;
        testCase.assertEqual(SDCD_R, SDCD_EXP)
      end

      % Size 0-->1
      ib_r      = [];
      ib_s      = [3];
      ib_r_asgn = [true];
      ib_exp    = [3];
      test()
      ib_r      = [];
      ib_s      = [3];
      ib_r_asgn = [1];
      ib_exp    = [3];
      test()

      % Overwrite 1 of 3.
      ib_r      = [];
      ib_s      = [3];
      ib_r_asgn = [true];
      ib_exp    = [3];
      test()

      % Overwrite 1 of 3 preexisting, add 1.
      ib_r      = [1,2,3];
      ib_s      = [4,5];
      ib_r_asgn = [2,4];
      ib_exp    = [1,4,3,5];
      test()
      ib_r      = [1,2,3];
      ib_s      = [4,5];
      ib_r_asgn = logical([0,1,0,1]);
      ib_exp    = [1,4,3,5];
      test()
    end



    function test_plus(testCase)
      SDCD_1     = bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1));
      SDCD_2     = bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1));
      testCase.assertEqual(...
        SDCD_1 + SDCD_2, ...
        bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1)))

      SDCD_1     = bicas.proc.L1L2.SdChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      SDCD_2     = bicas.proc.L1L2.SdChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      testCase.assertEqual(...
        SDCD_1 + SDCD_2, ...
        bicas.proc.L1L2.SdChannelData([2,5; 8,11; 14,NaN], logical([1; 1; 0])))
    end



    function test_minus(testCase)
      SDCD_1     = bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1));
      SDCD_2     = bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1));
      testCase.assertEqual(...
        SDCD_1 - SDCD_2, ...
        bicas.proc.L1L2.SdChannelData(zeros(0, 2), false(0, 1)))

      SDCD_1     = bicas.proc.L1L2.SdChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      SDCD_2     = bicas.proc.L1L2.SdChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      testCase.assertEqual(...
        SDCD_1 - SDCD_2, ...
        bicas.proc.L1L2.SdChannelData([0,1; 2,3; 4,NaN], logical([1; 1; 0])))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Get test SDCD with test data which is indexed before being used to create
    % SDCD.
    function Sdcd = get_SDCD(testCase, ib)
      SAMPLES_AR = [1,2; 3,4; 5,6; 7,8; 9,10];
      VSQB_AR    = logical([0;1;0;1;0]);

      Sdcd = bicas.proc.L1L2.SdChannelData(SAMPLES_AR(ib, :), VSQB_AR(ib, :));
    end



  end    % methods(Access=private)



end
