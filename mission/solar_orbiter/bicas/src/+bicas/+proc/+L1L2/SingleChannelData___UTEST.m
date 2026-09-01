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



    function test_constructor_size_bWholeRowIsNan(T)
      Schd = bicas.proc.L1L2.SingleChannelData(zeros(0, 0), false(0, 1));
      T.assertEqual(size(Schd),  [0, 1])
      T.assertEqual(Schd.bWholeRowIsNan, false(0, 1))

      Schd = bicas.proc.L1L2.SingleChannelData(zeros(0, 3), false(0, 1));
      T.assertEqual(size(Schd),  [0, 1])
      T.assertEqual(Schd.bWholeRowIsNan, false(0, 1))

      Schd = bicas.proc.L1L2.SingleChannelData(zeros(1, 3), false(1, 1));
      T.assertEqual(size(Schd),  [1, 1])
      T.assertEqual(Schd.bWholeRowIsNan, false(1, 1))



      SAMPLES_AR = [1,2; 3,4; 5,NaN; NaN,8; NaN,NaN];
      VSIB_AR    = logical([0; 1; 0; 1; 0]);

      Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR, VSIB_AR);
      T.assertEqual(Schd.samplesAr,      SAMPLES_AR)
      T.assertEqual(Schd.vsibAr,         VSIB_AR)
      T.assertEqual(size(Schd),          [5, 1])
      T.assertEqual(Schd.bWholeRowIsNan, logical([0; 0; 0; 0; 1]))
    end



    function test_subsref(T)
      SCHD_123 = T.get_SCHD([1,2,3]);
      SCHD_13  = T.get_SCHD([1,3]);
      SCHD_31  = T.get_SCHD([3,1]);
      SCHD_2   = T.get_SCHD([2]);
      SCHD_    = T.get_SCHD([]);

      T.assertEqual(SCHD_123, SCHD_123(1:3))
      T.assertEqual(SCHD_123, SCHD_123(logical([1,1,1])))

      T.assertEqual(SCHD_,   SCHD_([]))
      T.assertEqual(SCHD_,   SCHD_(logical([])))

      T.assertEqual(SCHD_13, SCHD_123(logical([1,0,1])))
      T.assertEqual(SCHD_13, SCHD_123([1,3]))

      T.assertEqual(SCHD_31, SCHD_123([3,1]))

      T.assertEqual(SCHD_2,  SCHD_123(logical([0,1,0])))
      T.assertEqual(SCHD_2,  SCHD_123([2]))

      T.assertEqual(SCHD_,   SCHD_123(logical([0,0,0])))
      T.assertEqual(SCHD_,   SCHD_123([]))
    end



    function test_subsasgn(T)
      function test()
        % R = Receiver?!
        % S = Sender?!
        SCHD_R            = T.get_SCHD(ib_r);
        SCHD_S            = T.get_SCHD(ib_s);
        SCHD_EXP          = T.get_SCHD(ib_exp);

        % Test subsasgn().
        % NOTE: SCHD is a value class, but subsasgn() still modifies the object
        %       in-place.
        SCHD_R(ib_r_asgn) = SCHD_S;

        T.assertEqual(SCHD_R, SCHD_EXP)
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



    function test_plus(T)
      SCHD_1     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      T.assertEqual(...
        SCHD_1 + SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1)))

      SCHD_1     = bicas.proc.L1L2.SingleChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      T.assertEqual(...
        SCHD_1 + SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData([2,5; 8,11; 14,NaN], logical([1; 1; 0])))
    end



    function test_minus(T)
      SCHD_1     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1));
      T.assertEqual(...
        SCHD_1 - SCHD_2, ...
        bicas.proc.L1L2.SingleChannelData(zeros(0, 2), false(0, 1)))

      SCHD_1     = bicas.proc.L1L2.SingleChannelData([1,3; 5,7; 9,NaN], logical([0; 1; 0]));
      SCHD_2     = bicas.proc.L1L2.SingleChannelData([1,2; 3,4; 5,6  ], logical([1; 0; 0]));
      T.assertEqual(...
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
    function Schd = get_SCHD(T, ib)
      SAMPLES_AR = [1,2; 3,4; 5,6; 7,8; 9,10];
      VSIB_AR    = logical([0;1;0;1;0]);

      Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR(ib, :), VSIB_AR(ib, :));
    end



  end    % methods(Access=private)



end
