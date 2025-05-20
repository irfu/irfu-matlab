%
% matlab.unittest automatic test code for bicas.proc.L1L2.SdChannelsData.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SingleChannelData___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor_size_bWholeRowIsNan(testCase)
      Schd = bicas.proc.L1L2.SingleChannelData(zeros(0, 0), false(0, 1));
      testCase.assertEqual(size(Schd),  [0, 1])
      testCase.assertEqual(Schd.bWholeRowIsNan, false(0, 1))

      Schd = bicas.proc.L1L2.SingleChannelData(zeros(0, 3), false(0, 1));
      testCase.assertEqual(size(Schd),  [0, 1])
      testCase.assertEqual(Schd.bWholeRowIsNan, false(0, 1))

      Schd = bicas.proc.L1L2.SingleChannelData(zeros(1, 3), false(1, 1));
      testCase.assertEqual(size(Schd),  [1, 1])
      testCase.assertEqual(Schd.bWholeRowIsNan, false(1, 1))



      SAMPLES_AR = [1,2; 3,4; 5,NaN; NaN,8; NaN,NaN];
      VSIB_AR    = logical([0; 1; 0; 1; 0]);

      Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR, VSIB_AR);
      testCase.assertEqual(Schd.samplesAr,      SAMPLES_AR)
      testCase.assertEqual(Schd.vsibAr,         VSIB_AR)
      testCase.assertEqual(size(Schd),          [5, 1])
      testCase.assertEqual(Schd.bWholeRowIsNan, logical([0; 0; 0; 0; 1]))
    end



    function test_subsref(testCase)
      SCHD_123 = testCase.get_SCHD([1,2,3]);
      SCHD_13  = testCase.get_SCHD([1,3]);
      SCHD_31  = testCase.get_SCHD([3,1]);
      SCHD_2   = testCase.get_SCHD([2]);
      SCHD_    = testCase.get_SCHD([]);

      testCase.assertEqual(SCHD_123, SCHD_123(1:3))
      testCase.assertEqual(SCHD_123, SCHD_123(logical([1,1,1])))

      testCase.assertEqual(SCHD_,   SCHD_([]))
      testCase.assertEqual(SCHD_,   SCHD_(logical([])))

      testCase.assertEqual(SCHD_13, SCHD_123(logical([1,0,1])))
      testCase.assertEqual(SCHD_13, SCHD_123([1,3]))

      testCase.assertEqual(SCHD_31, SCHD_123([3,1]))

      testCase.assertEqual(SCHD_2,  SCHD_123(logical([0,1,0])))
      testCase.assertEqual(SCHD_2,  SCHD_123([2]))

      testCase.assertEqual(SCHD_,   SCHD_123(logical([0,0,0])))
      testCase.assertEqual(SCHD_,   SCHD_123([]))
    end



    function test_subasgn(testCase)
      function test()
        % R = Receiver?!
        % S = Sender?!
        SCHD_R            = testCase.get_SCHD(ib_r);
        SCHD_S            = testCase.get_SCHD(ib_s);
        SCHD_EXP          = testCase.get_SCHD(ib_exp);
        SCHD_R(ib_r_asgn) = SCHD_S;
        testCase.assertEqual(SCHD_R, SCHD_EXP)
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
      SCHD_1     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      testCase.assertEqual(...
        SCHD_1 + SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1)))

      SCHD_1     = bicas.proc.L1L2.SingleChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      testCase.assertEqual(...
        SCHD_1 + SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData([2,5; 8,11; 14,NaN], logical([1; 1; 0])))
    end



    function test_minus(testCase)
      SCHD_1     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      testCase.assertEqual(...
        SCHD_1 - SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1)))

      SCHD_1     = bicas.proc.L1L2.SingleChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      testCase.assertEqual(...
        SCHD_1 - SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData([0,1; 2,3; 4,NaN], logical([1; 1; 0])))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Get test SCHD with test data which is indexed before being used to create
    % SCHD.
    function Schd = get_SCHD(testCase, ib)
      SAMPLES_AR = [1,2; 3,4; 5,6; 7,8; 9,10];
      VSIB_AR    = logical([0;1;0;1;0]);

      Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR(ib, :), VSIB_AR(ib, :));
    end



  end    % methods(Access=private)



end
