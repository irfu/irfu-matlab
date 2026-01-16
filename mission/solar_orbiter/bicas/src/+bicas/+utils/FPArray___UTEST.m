%
% matlab.unittest automatic test code for class bicas.utils.FPArray.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FPArray___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor(T)
      Fpa = bicas.utils.FPArray(      [1, 2, -1], 'FILL_VALUE', -1);
      T.assertEqual(Fpa.array(-2),    [1, 2, -2])
      T.assertEqual(Fpa.fpAr, logical([0, 0,  1]))

      Fpa = bicas.utils.FPArray(      [1, 2, -1], 'FILL_POSITIONS', [false, false, true]);
      T.assertEqual(Fpa.array(-2),    [1, 2, -2])
      T.assertEqual(Fpa.fpAr, logical([0, 0,  1]))

      % NO_FILL_POSITIONS
      Fpa = bicas.utils.FPArray(       [1, 2, -1]);
      T.assertEqual(Fpa.array(NaN),    [1, 2, -1])
      T.assertEqual(Fpa.fpAr,  logical([0, 0,  0]))
      Fpa = bicas.utils.FPArray(       [1, 2, -1], 'NO_FILL_POSITIONS');
      T.assertEqual(Fpa.array(NaN),    [1, 2, -1])
      T.assertEqual(Fpa.fpAr,  logical([0, 0,  0]))

      Fpa = bicas.utils.FPArray(      [ 1,  2, -1], 'ONLY_FILL_POSITIONS');
      T.assertEqual(Fpa.array(-2),    [-2, -2, -2])
      T.assertEqual(Fpa.fpAr, logical([ 1,  1,  1]))
    end



    function test_misc(T)
      % ===========================================================
      % Test legal & illegal read & write to *properties* (public &
      % private)
      % ===========================================================
      Fpa = bicas.utils.FPArray([]);

      % fpAr
      % NOTE: .fpAr is a READ-only property.
      function fpAr_assign_fail()
        Fpa.fpAr = logical([]);
      end
      T.verifyError(@() fpAr_assign_fail(), ?MException)
      T.verifyEqual(Fpa.fpAr, logical([]))

      % dataAr
      % NOTE: .dataAr is a private property.
      function dataAr_assign_fail()
        Fpa.dataAr = [];
      end
      T.verifyError(@() dataAr_assign_fail(), ?MException)
      if 0
        % Read dataAr (private property ==> should fail).
        % NOTE: Does not work for unknown reason!
        T.verifyError(@() (Fpa.dataAr), ?MException)
      end


      %=============
      % 0x0, double
      %=============
      if 1
        Fpa1 = bicas.utils.FPArray([]);
        Fpa2 = bicas.utils.FPArray([]);
        %
        T.test_equality(Fpa1, Fpa2, -99)
        T.verifyEqual(Fpa1.array(999), [])
        T.verifyEqual(Fpa1.fpAr,       logical([]))
      end

      %===========
      % 1D, int64
      %===========
      % Fill positions have different values.
      if 1
        Fpa1 = bicas.utils.FPArray(...
          int64([1, 2, 3, -8]), 'FILL_VALUE', int64(-8));
        Fpa2 = bicas.utils.FPArray(...
          int64([1, 2, 3, -9]), 'FILL_POSITIONS', logical([0, 0, 0, 1]));
        %
        T.test_equality(Fpa1, Fpa2, int64(-99))
        T.verifyEqual(...
          Fpa1.array(int64(-9)), int64([1, 2, 3, -9]))
        T.verifyEqual(...
          Fpa1.fpAr,             logical([0, 0, 0, 1]))
      end

      %==========
      % 2D, char
      %==========
      % Fill positions have different values.
      if 1
        Fpa1 = bicas.utils.FPArray(...
          ['abX'; 'dXf'], 'FILL_VALUE', 'X');
        Fpa2 = bicas.utils.FPArray(...
          ['abY'; 'dYf'], 'FILL_POSITIONS', logical([0, 0, 1; 0, 1, 0]));
        %
        T.test_equality(Fpa1, Fpa2, 'Z')
      end

      %============
      % 2D, double
      %============
      if 1
        fv1 = -1;
        fv2 = -2;
        Fpa1 = bicas.utils.FPArray(...
          [1, 2, fv1; 4, fv1, 6], 'FILL_VALUE', fv1);
        Fpa2 = bicas.utils.FPArray(...
          [1, 2,  -2; 4,  -5, 6], 'FILL_POSITIONS', ...
          logical([0, 0, 1; 0, 1, 0]));
        %
        T.test_equality(Fpa1, Fpa2, fv2)
        T.verifyEqual(...
          Fpa1.array(fv2), [1, 2,-2; 4,-2, 6])
        T.verifyEqual(...
          Fpa1.fpAr,       logical([0,0,1;0,1,0]))
      end

      %=====================================
      % Test NaN, +Inf, -Inf on minimal FPA
      % 1D, double
      % Both types of constructor calls.
      %=====================================
      if 1
        for x = [9, NaN, Inf, -Inf]               % Non-fill value.
          for fv1 = [-1, NaN, Inf, -Inf]        % Fill value when creating FPA.
            for fv2 = [-2, NaN, Inf, -Inf]    % Fill value when retrieving data from FPA.
              if ~isequaln(x, fv1)
                % fprintf('x=%g; fv1=%g; fv2=%g\n', x, fv1, fv2)

                Fpa1 = bicas.utils.FPArray(...
                  [x, fv1], 'FILL_VALUE', fv1);
                Fpa2 = bicas.utils.FPArray(...
                  [x,   x], 'FILL_POSITIONS', ...
                  logical([0, 1]));
                %
                T.test_equality(Fpa1, Fpa2, fv2)
                T.verifyEqual(...
                  Fpa1.array(fv2), ...
                  [x, fv2])
                T.verifyEqual(...
                  Fpa1.fpAr, ...
                  logical([0,1]))
              end
            end
          end
        end
      end

    end    % Function



    %         function test_set_FP(T)
    %             function test(ar, fv)
    %                 Fpa1  = T.Fpa(ar, fv);
    %                 fpAr0 = Fpa1.fpAr;
    %
    %                 Fpa1.set_FP();
    %                 fpAr1 = Fpa1.fpAr;
    %                 % Assert that fpAr has NOT changed.
    %                 T.assertEqual(fpAr0, fpAr1)
    %
    %                 Fpa2  = Fpa1.set_FP();
    %                 fpAr2 = Fpa2.fpAr;
    %                 % Assert that fpAr has changed/been set.
    %                 T.assertTrue(all(fpAr2, 'all'))
    %             end
    %
    %             test([], nan)
    %             test([1,2; 3,4], nan)
    %         end



    function test_convert(T)
      % Convert to type that forbids negative values, while input used
      % negative FV.
      Fpa1 = bicas.utils.FPArray([0,1,-1,2,3], 'FILL_VALUE', -1);
      Fpa2   = Fpa1.convert(@(x) (uint16(x)), 99);
      dataAr = Fpa2.array(uint16(999));
      T.verifyEqual(dataAr, uint16([0,1,999,2,3]))

      % Operation that raises error for NaN.
      % Convert to type that forbids NaN, while input used NaN as FV.
      Fpa1 = bicas.utils.FPArray([0, 3, NaN], 'FILL_VALUE', NaN);
      Fpa2   = Fpa1.convert(@(x) logical(~x), 99);
      dataAr = Fpa2.array(true);
      T.verifyEqual(dataAr, [true, false, true])
    end



    function test_cast(T)
      % double FPA --> logical FPA
      Fpa1   = bicas.utils.FPArray([0, 1, NaN], 'FILL_VALUE', NaN);
      Fpa2   = Fpa1.cast('logical');
      dataAr = Fpa2.array(true);
      T.verifyEqual(dataAr, [false, true, true])

      % logical FPA --> double FPA
      Fpa1   = bicas.utils.FPArray(...
        [false, true, false], 'FILL_POSITIONS', [false, false, true]);
      Fpa2   = Fpa1.cast('double');
      dataAr = Fpa2.array(NaN);
      T.verifyEqual(dataAr, [0, 1, NaN])

      % .cast() using automatic specification of FV (not specifying
      % argument).
      for mcCa1 = bicas.utils.FPArray.MC_NUMERIC_CA(:)'
        mc1 = mcCa1{1};
        for mcCa2 = bicas.utils.FPArray.MC_NUMERIC_CA(:)'
          mc2 = mcCa2{1};

          Fpa1 = bicas.utils.FPArray(...
            cast([false, false], mc1), 'FILL_POSITIONS', [false, true]);
          % Automatically set FV.
          Fpa2 = Fpa1.cast(mc2);
          % Set FV explicitly.
          Fpa3 = Fpa1.cast(mc2, cast(false, mc1));
        end
      end

    end



    function test_complement(T)
      Fpa1    = T.Fpa([], NaN);
      Fpa2    = T.Fpa([], NaN);
      ExpFpa3 = T.Fpa([],  -9);
      Fpa3    = Fpa1.complement(Fpa2);
      T.verifyTrue(Fpa3 == ExpFpa3);

      Fpa1    = T.Fpa([1, NaN;   3, NaN], NaN);
      Fpa2    = T.Fpa([1,   2; NaN, NaN], NaN);
      ExpFpa3 = T.Fpa([1,   2;   3,  -9],  -9);
      Fpa3 = Fpa1.complement(Fpa2);
      T.verifyTrue(Fpa3 == ExpFpa3);
    end



    function test_NFP_1D_array(T)
      % 0x0 --> 0x1
      Fpa1  = T.Fpa([], NaN);
      actAr = Fpa1.NFP_1D_array();
      T.assertEqual(actAr, zeros(0, 1))

      % Iterate over (1) different 1D vectors, in (2) different dimensions.
      for ca = {...
          {zeros(0,1),      zeros(0, 1)}, ...
          {[3, NaN, 4],     [3; 4]}, ...
          {[NaN, NaN], NaN, zeros(0, 1)} ...
          }'
        v1    = ca{1}{1};
        expAr = ca{1}{2};

        for iDim = 0:2
          v2 = permute(v1(:), wshift(1, 1:3, -iDim));

          Fpa1  = T.Fpa(v2, NaN);
          actAr = Fpa1.NFP_1D_array();
          T.assertEqual(actAr, expAr)
        end
      end

      % 2x3 --> 0x1
      Fpa1  = T.Fpa(NaN(2,3), NaN);
      actAr = Fpa1.NFP_1D_array();
      T.assertEqual(actAr, zeros(0,1))

      % 2x3 --> 3x1
      Fpa1  = T.Fpa([NaN, 2, NaN; 4, NaN, 6], NaN);
      actAr = Fpa1.NFP_1D_array();
      T.assertEqual(actAr, [4; 2; 6])
    end



    function test_ensure_NFP(T)
      Fpa1   = T.Fpa([], NaN);
      ExpFpa = T.Fpa([], NaN);
      ActFpa = Fpa1.ensure_NFP(-1);
      T.assertEqual(ActFpa, ExpFpa)

      Fpa1   = T.Fpa([3, NaN], NaN);
      ExpFpa = T.Fpa([3,  -1], NaN);
      ActFpa = Fpa1.ensure_NFP(-1);
      T.assertEqual(ActFpa, ExpFpa)
    end



    function test_subsref(T)
      % ===
      % 0x0
      % ===
      Fpa1 = T.Fpa([], NaN);
      Fpa2 = Fpa1(:);    % 0x0 --> 0x1
      T.assertEqual(size(Fpa2), [0, 1])

      Fpa1 = T.Fpa([], NaN);
      Fpa2 = Fpa1(:, :);    % 0x0 --> 0x0
      T.assertEqual(size(Fpa2), [0, 0])

      % ==
      % 2D
      % ==
      Fpa1 = T.Fpa([1,2; 3,4; 5,NaN], NaN);

      Fpa2 = Fpa1(3, 1:2);    % 2 values, one range.
      T.assertEqual(Fpa2, T.Fpa([5, NaN], NaN))

      % ===============
      % Linear indexing
      % ===============
      Fpa1 = T.Fpa([1,2; 3,4; 5,NaN], NaN);
      Fpa2 = Fpa1(:);
      Fpa3 = T.Fpa([1,3,5, 2,4,NaN]', NaN);
      T.assertEqual(Fpa2, Fpa3)

      Fpa1 = T.Fpa([1,2; 3,4; 5,NaN], NaN);
      Fpa2 = Fpa1([3,5,6]);
      Fpa3 = T.Fpa([5,4,NaN], NaN);
      T.assertEqual(Fpa2, Fpa3)

      % ================
      % Logical indexing
      % ================
      Fpa1 = T.Fpa([1,2; 3,4; 5,NaN], NaN);
      Fpa2 = Fpa1(logical([1,0; 0,0; 0,1]));    % 2D --> 1D
      Fpa3 = T.Fpa([1, NaN]', NaN);
      T.assertEqual(Fpa2, Fpa3)
    end



    function test_subsasgn_exact_size(T)
      % Copies FPAs, which requires FPA to be a non-handle class.
      assert(~isa(bicas.utils.FPArray.floatNan2logical([]), 'handle'))

      % Test assigning the wrong type: single to double. ==> Error
      Fpa1 = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa2 = T.Fpa(single([9]'), single(NaN));
      function test_assign_FPA_fail()
        Fpa1(1,1) = Fpa2;
      end
      T.verifyError(@() (test_assign_FPA_fail()), ?MException)
      function test_assign_array_fail()
        Fpa1(1,1) = Fpa2.array(NaN);
      end
      T.verifyError(@() (test_assign_array_fail()), ?MException)

      % Test keep both FV and non-FV.
      % Test overwrite FV/non-FV with FV/non-FV.
      Fpa1 = T.Fpa([1, NaN, 2, NaN,   3, NaN], NaN);
      Fpa2 = T.Fpa(        [8,   9, NaN, NaN], NaN);
      Fpa3 = T.Fpa([1, NaN, 8,   9, NaN, NaN], NaN);
      Fpa1(3:6) = Fpa2;
      T.assertEqual(Fpa1, Fpa3)

      % ================
      % Regular indexing
      % ================
      % Assign 0x0 to 0x0 over entire range.
      Fpa1a = T.Fpa([], NaN);
      Fpa1b = Fpa1a;

      Fpa2 = T.Fpa([], NaN);
      Fpa3 = T.Fpa([], NaN);
      Fpa1a(:, :, :, :) = Fpa2;
      T.assertEqual(Fpa1a, Fpa3)

      %             Fpa1b(:, :, :, :) = Fpa2.array(NaN);
      %             T.assertEqual(Fpa1b, Fpa3)

      % Assign single element.
      Fpa1 = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa2 = T.Fpa([9], NaN);
      Fpa3 = T.Fpa([1,2,3; 4,5,9], NaN);
      Fpa1(2, 3) = Fpa2;
      T.verifyEqual(Fpa1, Fpa3)

      % Assign entire row.
      Fpa1 = T.Fpa([1,2,NaN; 4,NaN,  6], NaN);
      Fpa2 = T.Fpa([9,8,NaN], NaN);
      Fpa3 = T.Fpa([1,2,NaN; 9,  8,NaN], NaN);
      Fpa1(2, :) = Fpa2;
      T.verifyEqual(Fpa1, Fpa3)

      % Assign entire dimensions.
      Fpa1 = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa2 = T.Fpa([9,8,7; 6,5,4], NaN);
      Fpa3 = Fpa2;
      Fpa1(:, :) = Fpa2;
      T.verifyEqual(Fpa1, Fpa3)

      % end
      % scalar input --> vector output
      Fpa1 = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa2 = T.Fpa([9], NaN);
      Fpa3 = T.Fpa([1,9,3; 4,9,6], NaN);
      Fpa1(1:end, end-1) = Fpa2;
      T.verifyEqual(Fpa1, Fpa3)

      % ===============
      % Linear indexing
      % ===============
      % Assign non-rectangular subset of elements in 2D FPA.
      Fpa1 = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa2 = T.Fpa([9,8,7]', NaN);
      Fpa3 = T.Fpa([9,2,7; 4,8,6], NaN);
      Fpa1([1,4,5]) = Fpa2;
      T.verifyEqual(Fpa1, Fpa3)

      % ================
      % Logical indexing
      % ================
      Fpa1a = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa1b = Fpa1a;
      Fpa2 = T.Fpa([9,8,7]', NaN);
      Fpa3 = T.Fpa([9,2,7; 4,8,6], NaN);
      Fpa1a(logical([1,0,1; 0,1,0])) = Fpa2;
      T.verifyEqual(Fpa1a, Fpa3)

      %             Fpa1b(logical([1,0,1; 0,1,0])) = Fpa2.array(NaN);
      %             T.verifyEqual(Fpa1b, Fpa3)

      % =======================
      % Logical linear indexing
      % =======================
      % Index = 1D vector of logical.
      Fpa1a = T.Fpa([1,2,3; 4,5,6], NaN);
      Fpa1b = Fpa1a;
      Fpa2 = T.Fpa([9,8,7]', NaN);
      Fpa3 = T.Fpa([9,2,7; 4,8,6], NaN);
      Fpa1a(logical([1,0,0, 1,1,0]')) = Fpa2;
      T.verifyEqual(Fpa1a, Fpa3)

      %             Fpa1b(logical([1,0,0, 1,1,0]')) = Fpa2.array(NaN);
      %             T.verifyEqual(Fpa1b, Fpa3)
    end



    function test_subsasgn_different_size(T)
      % Copies FPAs, which requires FPA to be a non-handle class.
      assert(~isa(bicas.utils.FPArray.floatNan2logical([]), 'handle'))

      % Assign scalar to matrix.
      Fpa1a = T.Fpa([1,NaN,3; 4,5,NaN], NaN);
      Fpa1b = Fpa1a;
      Fpa2  = T.Fpa([9], NaN);
      Fpa3  = T.Fpa([9,9,9; 9,9,9], NaN);
      Fpa1a(:, :) = Fpa2;
      T.verifyEqual(Fpa1a, Fpa3)
      %             Fpa1b(:, :) = Fpa2.array(NaN);
      %             T.verifyEqual(Fpa1b, Fpa3)

    end


    function test_size_ndims(T)
      % IMPLEMENTATION NOTE: In a sense, this does not only test the code,
      % but also the author's understanding of "overloading" with a method
      % size().

      TEST_DATA_CA = {...
        ones(0, 0), ...
        ones(1, 1), ...
        ones(1, 1, 0), ...
        ones(3, 1), ...      % 1D
        ones(1, 3), ...      % 1D, in second dimension.
        ones(4, 5), ...      % 2D
        ones(4, 5, 6), ...   % 3D
        };

      for ca = TEST_DATA_CA(:)'
        v = ca{1};
        Fpa = bicas.utils.FPArray(v, 'FILL_VALUE', NaN);

        T.verifyEqual(isscalar(Fpa), isscalar(v))
        T.verifyEqual(size(Fpa),     size(v)    )
        T.verifyEqual(ndims(Fpa),    ndims(v)   )
        T.verifyEqual(isempty(Fpa),  isempty(v) )

        for iDim = 1:3
          T.verifyEqual(size(Fpa, iDim), size(v, iDim) )
        end
      end    % for
    end



    function test_lt_gt_le_ge(T)
      % PROPOSAL: Combine with tests for == and min().
      %   PRO: Should have synergies. Similar tests.
      %   PRO: There are relationships between <, >, <=, >=, == that one could use.
      %
      % PROBLEM: There are many combinations to test. Would ideally like
      %          to test every combination of:
      %   input 1/2 elements: lower number, higher number, NaN, -Inf, +Inf, FP
      %   input 1/2 FPA     scalar and non-scalar
      %   input 2   non-FPA scalar and non-scalar
      %   operation: <, >, <=, >=
      %   PROPOSAL: Automatically derive expected results from using non-FPAs.
      %   PROPOSAL: Specify results for any combination of inputs for <.
      %             Then derive the remaining ones automatically.
      %
      % PROPOSAL: Test different MATLAB classes.

      % Define FPAs which are used as input for all operators. Should
      % contain every case of element-wise comparison.
      % NOTE: Expected FPs will be the same, regardless of operation.
      Fpa1 = T.Fpa([1,  3, 2,   2, NaN], NaN);
      Fpa2 = T.Fpa([3,  1, 2, NaN,   2], NaN);
      EXP_FP = [false, false, false,  true,  true];

      % ===
      %  <
      % ===
      % FPA non-scalar < FPA non-scalar
      ActFpa3 = Fpa1 < Fpa2;
      ExpFpa3 = bicas.utils.FPArray(...
        [true,  false, false, false, false], 'FILL_POSITIONS', EXP_FP);
      T.assertEqual(ActFpa3, ExpFpa3)
      % Non-scalar FPA < scalar FPA
      ActFpa3 = Fpa1 < T.Fpa(2, NaN);
      ExpFpa3 = bicas.utils.FPArray(...
        [true,  false, false, false, false], 'FILL_POSITIONS', ...
        [false, false, false, false, true ]);
      T.assertEqual(ActFpa3, ExpFpa3)
      % Scalar FPA < non-scalar FPA -- Tests specific bugfix!
      ActFpa3 = T.Fpa(2, NaN) < [1,2,3,NaN];
      ExpFpa3 = bicas.utils.FPArray(...
        [false, false, true,  false], 'FILL_POSITIONS', ...
        [false, false, false, false]);
      T.assertEqual(ActFpa3, ExpFpa3)

      % ===
      %  >
      % ===
      % FPA non-scalar > FPA non-scalar
      ActFpa3 = Fpa1 > Fpa2;
      ExpFpa3 = bicas.utils.FPArray(...
        [false,  true, false, false, false], 'FILL_POSITIONS', EXP_FP);
      T.assertEqual(ActFpa3, ExpFpa3)
      % Non-scalar FPA > scalar FPA
      ActFpa3 = T.Fpa(2, NaN) > Fpa2;
      ExpFpa3 = bicas.utils.FPArray(...
        [false,  true, false, false, false], 'FILL_POSITIONS', ...
        [false, false, false, true,  false]);
      T.assertEqual(ActFpa3, ExpFpa3)

      % ====
      %  <=
      % ====
      ActFpa3 = Fpa1 <= Fpa2;
      ExpFpa3 = bicas.utils.FPArray(...
        [true,  false, true, false, false], 'FILL_POSITIONS', EXP_FP);
      T.assertEqual(ActFpa3, ExpFpa3)

      % ====
      %  >=
      % ====
      ActFpa3 = Fpa1 >= Fpa2;
      ExpFpa3 = bicas.utils.FPArray(...
        [false,  true, true, false, false], 'FILL_POSITIONS', EXP_FP);
      T.assertEqual(ActFpa3, ExpFpa3)
    end



    function test_plus(T)
      function test(Fpa1, Fpa2, ExpFpa)
        ActFpa12  = Fpa1 + Fpa2;
        ActFpa21  = Fpa2 + Fpa1;
        ActFpa12b = Fpa1 + Fpa2.array();
        ActFpa21b = Fpa2 + Fpa1.array();

        T.assertEqual(ActFpa12,  ExpFpa)
        T.assertEqual(ActFpa21,  ExpFpa)
        T.assertEqual(ActFpa12b, ExpFpa)
        T.assertEqual(ActFpa21b, ExpFpa)
      end

      % Non-scalar + non-scalar
      test(T.Fpa([2, 3], NaN), T.Fpa([4, 5], NaN), T.Fpa([6, 8], NaN))
      % Scalar + non-scalar
      test(T.Fpa([2],    NaN), T.Fpa([4; 5], NaN), T.Fpa([6; 7], NaN))
      % Scalar + scalar
      test(T.Fpa([2],    NaN), T.Fpa([4],    NaN), T.Fpa([6],    NaN))
    end



    function test_times(T)
      function test(Fpa1, Fpa2, ExpFpa)
        ActFpa12  = Fpa1 .* Fpa2;
        ActFpa21  = Fpa2 .* Fpa1;
        ActFpa12b = Fpa1 .* Fpa2.array();
        ActFpa21b = Fpa2 .* Fpa1.array();

        T.assertEqual(ActFpa12,  ExpFpa)
        T.assertEqual(ActFpa21,  ExpFpa)
        T.assertEqual(ActFpa12b, ExpFpa)
        T.assertEqual(ActFpa21b, ExpFpa)
      end

      % Non-scalar .* non-scalar
      test(T.Fpa([2, 3], NaN), T.Fpa([3, 4], NaN), T.Fpa([6, 12], NaN))
      % Scalar .* non-scalar
      test(T.Fpa([2],    NaN), T.Fpa([3; 4], NaN), T.Fpa([6;  8], NaN))
      % Scalar .* scalar
      test(T.Fpa([2],    NaN), T.Fpa([3],    NaN), T.Fpa([6],     NaN))
    end



    %         function test_or(T)
    %             function Fpa = Fpa(dataAr, fpAr)
    %                 Fpa = bicas.utils.FPArray(logical(dataAr), ...
    %                     'FILL_POSITIONS', logical(fpAr));
    %             end
    %
    %             function test(Fpa1, Fpa2, ExpFpa)
    %                 ActFpa12  = Fpa1 | Fpa2;
    %                 ActFpa21  = Fpa2 | Fpa1;
    %                 ActFpa12b = Fpa1 | Fpa2.array();
    %                 ActFpa21b = Fpa2 | Fpa1.array();
    %
    %                 T.assertEqual(ActFpa12,  ExpFpa)
    %                 T.assertEqual(ActFpa21,  ExpFpa)
    %                 T.assertEqual(ActFpa12b, ExpFpa)
    %                 T.assertEqual(ActFpa21b, ExpFpa)
    %             end
    %
    %             % Non-scalar | non-scalar
    %             test(...
    %                 Fpa([0 1 0 1; 0 1 0 1], [0 0 0 0; 1 1 1 1]), ...
    %                 Fpa([0 0 1 1; 0 0 1 1], [0 0 0 0; 1 1 1 1]), ...
    %                 Fpa([0 1 1 1; 1 0 0 0], [0 0 0 0; 1 1 1 1]))
    %
    %             % Scalar | non-scalar
    %             test(...
    %                 Fpa([0],    [0]), ...
    %                 Fpa([0; 1], [0; 0]), ...
    %                 Fpa([0; 1], [0; 0]))
    %
    %             % Scalar | scalar
    %             test(Fpa([1], [0]), Fpa([1], [0]), Fpa([1], [0]))
    %             test(Fpa([0], [0]), Fpa([0], [0]), Fpa([0], [0]))
    %         end



    function test_floatNan2logical(T)
      function test_element_illegal_fail(mc2)
        Fpa = bicas.utils.FPArray.floatNan2logical(...
          cast([2], mc2));
      end

      for mcCa = {'single', 'double'}'
        mc = mcCa{1};

        ExpFpa = bicas.utils.FPArray(...
          [false, true,  false], 'FILL_POSITIONS', ...
          [false, false,  true]);
        ActFpa = bicas.utils.FPArray.floatNan2logical(...
          cast([0, 1, NaN], mc));
        T.verifyTrue(ExpFpa == ActFpa)

        T.verifyError(@() test_element_illegal_fail(mc), ?MException)
      end
    end



    function test_floatNan2int(T)
      for floatMcCa = {'single', 'double'}'
        floatMc = floatMcCa{1};
        for intMcCa = {'uint8', 'int8'}'
          intMc = intMcCa{1};

          ExpFpa = bicas.utils.FPArray(...
            cast([    0,     1,     0], intMc), 'FILL_POSITIONS', ...
            [false, false,  true]);

          ActFpa = bicas.utils.FPArray.floatNan2int(...
            cast([0, 1, NaN], floatMc), intMc);

          T.verifyTrue(ExpFpa == ActFpa)
        end
      end
    end



    function test_cat_vertcat_horzcat(T)
      Fpa1 = T.Fpa([], nan);
      ActFpa = [Fpa1];
      ExpFpa = Fpa1;
      T.assertEqual(ActFpa, ExpFpa)

      Fpa1 = T.Fpa([1, 2], nan);
      Fpa2 = T.Fpa([3, 4], nan);

      ActFpa = [Fpa1, Fpa2];
      ExpFpa = T.Fpa([1, 2, 3, 4], nan);
      T.assertEqual(ActFpa, ExpFpa)

      ActFpa = [Fpa1; Fpa2];
      ExpFpa = T.Fpa([1, 2; 3, 4], nan);
      T.assertEqual(ActFpa, ExpFpa)

      ActFpa = cat(1, Fpa1, Fpa2);
      ExpFpa = T.Fpa([1, 2; 3, 4], nan);
      T.assertEqual(ActFpa, ExpFpa)

      ActFpa = cat(2, Fpa1, Fpa2);
      ExpFpa = T.Fpa([1, 2, 3, 4], nan);
      T.assertEqual(ActFpa, ExpFpa)
    end



    function test_repmat(T)
      EMPTY_FPA = T.Fpa([], NaN);

      ActFpa  = repmat(EMPTY_FPA, [2,3]);
      T.assertEqual(ActFpa, EMPTY_FPA)
      ActFpa  = repmat(EMPTY_FPA, 2, 3);
      T.assertEqual(ActFpa, EMPTY_FPA)

      DATA_AR = [1 2 3; 4 5 NaN];

      FPA_0  = T.Fpa(DATA_AR, NaN);
      ActFpa = repmat(FPA_0, 1);
      T.assertEqual(ActFpa, FPA_0)

      for ca = {{1}, {2}, {[2,3]}, {2,3}, {[2,3,4]}}
        repmatArgsCa = ca{1};

        T.assertEqual(...
          repmat(FPA_0, repmatArgsCa{:}), ...
          T.Fpa(repmat(DATA_AR, repmatArgsCa{:}), NaN) ...
          )
      end
    end



    function test_min(T)
      Fpa1   = T.Fpa([1, 2, NaN, 4  ], NaN);
      Fpa2   = T.Fpa([2, 1, 3,   NaN], NaN);
      ExpFpa = T.Fpa([1, 1, NaN, NaN], NaN);
      ActFpa = Fpa1.min(Fpa2);
      T.assertEqual(ActFpa, ExpFpa)

      Fpa2    = T.Fpa([2], NaN);
      ExpFpa  = T.Fpa([1, 2, NaN, 2], NaN);
      ActFpa1 = Fpa1.min(Fpa2);
      ActFpa2 = Fpa2.min(Fpa1);
      T.assertEqual(ActFpa1, ExpFpa)
      T.assertEqual(ActFpa2, ExpFpa)
    end



    % NOTE: Only tests the method indirectly, and only by checking if code
    % does not crash.
    function test_getPropertyGroups(T)
      T.Fpa(ones(0,0), NaN);
      T.Fpa(ones(0,1), NaN);
      T.Fpa(ones(1,0), NaN);
      T.Fpa(ones(2,3,4), NaN);
      T.Fpa([NaN,inf,-inf], -1);
      T.Fpa([0,1,2; 3,4,5], NaN);
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Do multiple tests for FPAs which should be EQUAL.
    function test_equality(T, fpa1, fpa2, fv)
      % NOTE: For tests to be truly meaningful, implementation-dependent
      % fill position values should different.

      % Ensure that method "eq" is called.
      r = (fpa1 == fpa2);
      T.assertTrue(isscalar(r))
      T.verifyTrue(r)

      % Ensure that method "ne" is called.
      r = (fpa1 ~= fpa2);
      T.assertTrue(isscalar(r))
      T.verifyFalse(r)

      % Ensure that method "isequaln" is called.
      r = isequaln(fpa1, fpa2);
      T.assertTrue(isscalar(r))
      T.verifyTrue(r)

      % Verify equality through the methods for obtaining data.
      T.verifyEqual(...
        fpa1.array(fv), ...
        fpa2.array(fv))
      T.verifyEqual(...
        fpa1.fpAr, ...
        fpa2.fpAr)
    end



  end    % methods(Access=private)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function Fpa = Fpa(dataAr, fv)
      Fpa = bicas.utils.FPArray(dataAr, 'FILL_VALUE', fv);
    end



  end    % methods(Static, Access=private)



end
