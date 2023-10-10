%
% matlab.unittest automatic test code for class bicas.utils.FillPositionsArray.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FillPositionsArray___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)
        
        
        
        function test_constructor(testCase)
            Fpa = bicas.utils.FillPositionsArray(  [1, 2, -1], 'FILL_VALUE', -1);
            testCase.assertEqual(Fpa.get_data(-2), [1, 2, -2])
            testCase.assertEqual(Fpa.fpAr, logical([0, 0,  1]))

            Fpa = bicas.utils.FillPositionsArray(  [1, 2, -1], 'FILL_POSITIONS', [false, false, true]);
            testCase.assertEqual(Fpa.get_data(-2), [1, 2, -2])
            testCase.assertEqual(Fpa.fpAr, logical([0, 0,  1]))

            Fpa = bicas.utils.FillPositionsArray(   [1, 2, -1], 'NO_FILL_POSITIONS');
            testCase.assertEqual(Fpa.get_data(NaN), [1, 2, -1])
            testCase.assertEqual(Fpa.fpAr,  logical([0, 0,  0]))

            Fpa = bicas.utils.FillPositionsArray(  [1,   2, -1], 'ONLY_FILL_POSITIONS');
            testCase.assertEqual(Fpa.get_data(-2), [-2, -2, -2])
            testCase.assertEqual(Fpa.fpAr, logical([ 1,  1,  1]))
        end



        function test_misc(testCase)
            import bicas.utils.FillPositionsArray___UTEST.test_equality

            % ===========================================================
            % Test legal & illegal read & write to *properties* (public &
            % private)
            % ===========================================================
            Fpa = bicas.utils.FillPositionsArray([], 'NO_FILL_POSITIONS');
            
            % fpAr
            % NOTE: .fpAr is a READ-only property.
            function fpAr_assign_fail()
                Fpa.fpAr = logical([]);
            end
            testCase.verifyError(@() fpAr_assign_fail(), ?MException)
            testCase.verifyEqual(Fpa.fpAr, logical([]))
            
            % dataAr
            % NOTE: .dataAr is a private property.
            function dataAr_assign_fail()
                Fpa.dataAr = [];
            end
            testCase.verifyError(@() dataAr_assign_fail(), ?MException)
            if 0
                % Read dataAr (should fail).
                % NOTE: Does not work, since class does not work.
                testCase.verifyError(@() (Fpa.dataAr), ?MException)
            end
            
            
            %=============
            % 0x0, double
            %=============
            if 1
                Fpa1 = bicas.utils.FillPositionsArray([], 'NO_FILL_POSITIONS');
                Fpa2 = bicas.utils.FillPositionsArray([], 'NO_FILL_POSITIONS');
                %
                test_equality(testCase, Fpa1, Fpa2, -99)
                testCase.verifyEqual(Fpa1.get_data(999), [])
                testCase.verifyEqual(Fpa1.fpAr,          logical([]))
            end

            %===========
            % 1D, int64
            %===========
            % Fill positions have different values.
            if 1
                Fpa1 = bicas.utils.FillPositionsArray(...
                    int64([1, 2, 3, -8]), 'FILL_VALUE', int64(-8));
                Fpa2 = bicas.utils.FillPositionsArray(...
                    int64([1, 2, 3, -9]), 'FILL_POSITIONS', logical([0, 0, 0, 1]));
                %
                test_equality(testCase, Fpa1, Fpa2, int64(-99))
                testCase.verifyEqual(...
                    Fpa1.get_data(int64(-9)), int64([1, 2, 3, -9]))
                testCase.verifyEqual(...
                    Fpa1.fpAr,                logical([0,0,0,1]))
            end
            
            %==========
            % 2D, char
            %==========
            % Fill positions have different values.
            if 1
                Fpa1 = bicas.utils.FillPositionsArray(...
                    ['abX'; 'dXf'], 'FILL_VALUE', 'X');
                Fpa2 = bicas.utils.FillPositionsArray(...
                    ['abY'; 'dYf'], 'FILL_POSITIONS', logical([0, 0, 1; 0, 1, 0]));
                %
                test_equality(testCase, Fpa1, Fpa2, 'Z')
            end
            
            %============
            % 2D, double
            %============
            if 1
                fv1 = -1;
                fv2 = -2;
                Fpa1 = bicas.utils.FillPositionsArray(...
                    [1, 2, fv1; 4, fv1, 6], 'FILL_VALUE', fv1);
                Fpa2 = bicas.utils.FillPositionsArray(...
                    [1, 2,  -2; 4,  -5, 6], 'FILL_POSITIONS', ...
                    logical([0, 0, 1; 0, 1, 0]));
                %
                test_equality(testCase, Fpa1, Fpa2, fv2)
                testCase.verifyEqual(...
                    Fpa1.get_data(fv2), [1, 2,-2; 4,-2, 6])
                testCase.verifyEqual(...
                    Fpa1.fpAr,          logical([0,0,1;0,1,0]))
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

                                Fpa1 = bicas.utils.FillPositionsArray(...
                                    [x, fv1], 'FILL_VALUE', fv1);
                                Fpa2 = bicas.utils.FillPositionsArray(...
                                    [x,   x], 'FILL_POSITIONS', ...
                                    logical([0, 1]));
                                %
                                test_equality(testCase, Fpa1, Fpa2, fv2)
                                testCase.verifyEqual(...
                                    Fpa1.get_data(fv2), ...
                                    [x, fv2])
                                testCase.verifyEqual(...
                                    Fpa1.fpAr, ...
                                    logical([0,1]))
                            end
                        end
                    end
                end
            end

        end    % Function
        
        
        
%         function test_set_FP(testCase)
%             import bicas.utils.FillPositionsArray___UTEST.Fpa
%             
%             function test(ar, fv)
%                 Fpa1  = Fpa(ar, fv);
%                 fpAr0 = Fpa1.fpAr;
%                 
%                 Fpa1.set_FP();
%                 fpAr1 = Fpa1.fpAr;
%                 % Assert that fpAr has NOT changed.
%                 testCase.assertEqual(fpAr0, fpAr1)
%                 
%                 Fpa2  = Fpa1.set_FP();
%                 fpAr2 = Fpa2.fpAr;
%                 % Assert that fpAr has changed/been set.
%                 testCase.assertTrue(all(fpAr2, 'all'))
%             end
%             
%             test([], nan)
%             test([1,2; 3,4], nan)
%         end



        function test_convert(testCase)
            % Convert to type that forbids negative values, while input used
            % negative FV.
            Fpa1 = bicas.utils.FillPositionsArray([0,1,-1,2,3], 'FILL_VALUE', -1);
            Fpa2   = Fpa1.convert(@(x) (uint16(x)), 99);
            dataAr = Fpa2.get_data(uint16(999));
            testCase.verifyEqual(dataAr, uint16([0,1,999,2,3]))
            
            % Operation that raises error for NaN.
            % Convert to type that forbids NaN, while input used NaN as FV.
            Fpa1 = bicas.utils.FillPositionsArray([0, 3, NaN], 'FILL_VALUE', NaN);
            Fpa2   = Fpa1.convert(@(x) logical(~x), 99);
            dataAr = Fpa2.get_data(true);
            testCase.verifyEqual(dataAr, [true, false, true])
        end
        
        
        
        function test_cast(testCase)
            % double FPA --> logical FPA
            Fpa1   = bicas.utils.FillPositionsArray([0, 1, NaN], 'FILL_VALUE', NaN);
            Fpa2   = Fpa1.cast('logical');
            dataAr = Fpa2.get_data(true);
            testCase.verifyEqual(dataAr, [false, true, true])
            
            % logical FPA --> double FPA
            Fpa1   = bicas.utils.FillPositionsArray(...
                [false, true, false], 'FILL_POSITIONS', [false, false, true]);
            Fpa2   = Fpa1.cast('double');
            dataAr = Fpa2.get_data(NaN);
            testCase.verifyEqual(dataAr, [0, 1, NaN])
            
            % .cast() using automatic specification of FV (not specifying
            % argument).
            for mcCa1 = bicas.utils.FillPositionsArray.MC_NUMERIC_CA(:)'
                mc1 = mcCa1{1};
                for mcCa2 = bicas.utils.FillPositionsArray.MC_NUMERIC_CA(:)'
                    mc2 = mcCa2{1};
                    
                    Fpa1 = bicas.utils.FillPositionsArray(...
                        cast([false, false], mc1), 'FILL_POSITIONS', [false, true]);
                    % Automatically set FV.
                    Fpa2 = Fpa1.cast(mc2);
                    % Set FV explicitly.
                    Fpa3 = Fpa1.cast(mc2, cast(false, mc1));
                end
            end
                
        end
        
        
        
        function test_complement(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            Fpa1    = Fpa([], NaN);
            Fpa2    = Fpa([], NaN);            
            ExpFpa3 = Fpa([],  -9);
            Fpa3    = Fpa1.complement(Fpa2);
            testCase.verifyTrue(Fpa3 == ExpFpa3);

            Fpa1    = Fpa([1, NaN;   3, NaN], NaN);
            Fpa2    = Fpa([1,   2; NaN, NaN], NaN);
            ExpFpa3 = Fpa([1,   2;   3,  -9],  -9);
            Fpa3 = Fpa1.complement(Fpa2);
            testCase.verifyTrue(Fpa3 == ExpFpa3);
        end
        
        
        
        function test_get_non_FP_data(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            % 0x0 --> 0x1
            Fpa1  = Fpa([], NaN);
            actAr = Fpa1.get_non_FP_data();
            testCase.assertEqual(actAr, zeros(0, 1))
            
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
%                     disp('---------')
%                     size(v1)
%                     size(v2)                    

                    Fpa1  = Fpa(v2, NaN);
                    actAr = Fpa1.get_non_FP_data();
                    testCase.assertEqual(actAr, expAr)
                end
            end
            
            % 2x3 --> 0x1
            Fpa1  = Fpa(NaN(2,3), NaN);
            actAr = Fpa1.get_non_FP_data();
            testCase.assertEqual(actAr, zeros(0,1))

            % 2x3 --> 3x1
            Fpa1  = Fpa([NaN, 2, NaN; 4, NaN, 6], NaN);
            actAr = Fpa1.get_non_FP_data();
            testCase.assertEqual(actAr, [4; 2; 6])
        end
        
        
        
        function test_subsref(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            % ===
            % 0x0
            % ===
            Fpa1 = Fpa([], NaN);
            Fpa2 = Fpa1(:);    % 0x0 --> 0x1
            testCase.assertEqual(size(Fpa2), [0, 1])
            
            Fpa1 = Fpa([], NaN);
            Fpa2 = Fpa1(:, :);    % 0x0 --> 0x0
            testCase.assertEqual(size(Fpa2), [0, 0])

            % ==
            % 2D
            % ==
            Fpa1 = Fpa([1,2; 3,4; 5,NaN], NaN);
            
            Fpa2 = Fpa1(3, 1:2);    % 2 values, one range.
            testCase.assertEqual(Fpa2, Fpa([5, NaN], NaN))
            
            % ===============
            % Linear indexing
            % ===============
            Fpa1 = Fpa([1,2; 3,4; 5,NaN], NaN);
            Fpa2 = Fpa1(:);
            Fpa3 = Fpa([1,3,5, 2,4,NaN]', NaN);
            testCase.assertEqual(Fpa2, Fpa3)
            
            Fpa1 = Fpa([1,2; 3,4; 5,NaN], NaN);
            Fpa2 = Fpa1([3,5,6]);
            Fpa3 = Fpa([5,4,NaN], NaN);
            testCase.assertEqual(Fpa2, Fpa3)
            
            % ================
            % Logical indexing
            % ================
            Fpa1 = Fpa([1,2; 3,4; 5,NaN], NaN);
            Fpa2 = Fpa1(logical([1,0; 0,0; 0,1]));    % 2D --> 1D
            Fpa3 = Fpa([1, NaN]', NaN);
            testCase.assertEqual(Fpa2, Fpa3)
        end
        
        
        
        function test_subsasgn_exact_size(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            % Copies FPAs, which requires FPA to be a non-handle class.
            assert(~isa(bicas.utils.FillPositionsArray.floatNan2logical([]), 'handle'))
            
            % Test assigning the wrong type: single to double. ==> Error
            Fpa1 = Fpa([1,2,3; 4,5,6], NaN);
            Fpa2 = Fpa(single([9]'), single(NaN));
            function test_assign_FPA_fail()
                Fpa1(1,1) = Fpa2;
            end
            testCase.verifyError(@() (test_assign_FPA_fail()), ?MException)            
            function test_assign_array_fail()
                Fpa1(1,1) = Fpa2.get_data(NaN);
            end
            testCase.verifyError(@() (test_assign_array_fail()), ?MException)            

            % Test keep both FV and non-FV.
            % Test overwrite FV/non-FV with FV/non-FV.
            Fpa1 = Fpa([1, NaN, 2, NaN,   3, NaN], NaN);
            Fpa2 = Fpa(        [8,   9, NaN, NaN], NaN);
            Fpa3 = Fpa([1, NaN, 8,   9, NaN, NaN], NaN);
            Fpa1(3:6) = Fpa2;
            testCase.assertEqual(Fpa1, Fpa3)

            % ================
            % Regular indexing
            % ================
            % Assign 0x0 to 0x0 over entire range.
            Fpa1a = Fpa([], NaN);
            Fpa1b = Fpa1a;
            
            Fpa2 = Fpa([], NaN);
            Fpa3 = Fpa([], NaN);
            Fpa1a(:, :, :, :) = Fpa2;
            testCase.assertEqual(Fpa1a, Fpa3)

%             Fpa1b(:, :, :, :) = Fpa2.get_data(NaN);
%             testCase.assertEqual(Fpa1b, Fpa3)

            % Assign single element.
            Fpa1 = Fpa([1,2,3; 4,5,6], NaN);
            Fpa2 = Fpa([9], NaN);
            Fpa3 = Fpa([1,2,3; 4,5,9], NaN);            
            Fpa1(2, 3) = Fpa2;
            testCase.verifyEqual(Fpa1, Fpa3)
            
            % Assign entire row.
            Fpa1 = Fpa([1,2,NaN; 4,NaN,  6], NaN);
            Fpa2 = Fpa([9,8,NaN], NaN);
            Fpa3 = Fpa([1,2,NaN; 9,  8,NaN], NaN);
            Fpa1(2, :) = Fpa2;
            testCase.verifyEqual(Fpa1, Fpa3)

            % Assign entire dimensions.
            Fpa1 = Fpa([1,2,3; 4,5,6], NaN);
            Fpa2 = Fpa([9,8,7; 6,5,4], NaN);
            Fpa3 = Fpa2;
            Fpa1(:, :) = Fpa2;
            testCase.verifyEqual(Fpa1, Fpa3)
            
            % end
            % scalar input --> vector output
            Fpa1 = Fpa([1,2,3; 4,5,6], NaN);
            Fpa2 = Fpa([9], NaN);
            Fpa3 = Fpa([1,9,3; 4,9,6], NaN);
            Fpa1(1:end, end-1) = Fpa2;
            testCase.verifyEqual(Fpa1, Fpa3)
            
            % ===============
            % Linear indexing
            % ===============
            % Assign non-rectangular subset of elements in 2D FPA.
            Fpa1 = Fpa([1,2,3; 4,5,6], NaN);
            Fpa2 = Fpa([9,8,7]', NaN);
            Fpa3 = Fpa([9,2,7; 4,8,6], NaN);
            Fpa1([1,4,5]) = Fpa2;
            testCase.verifyEqual(Fpa1, Fpa3)
            
            % ================
            % Logical indexing
            % ================
            Fpa1a = Fpa([1,2,3; 4,5,6], NaN);
            Fpa1b = Fpa1a;
            Fpa2 = Fpa([9,8,7]', NaN);
            Fpa3 = Fpa([9,2,7; 4,8,6], NaN);
            Fpa1a(logical([1,0,1; 0,1,0])) = Fpa2;
            testCase.verifyEqual(Fpa1a, Fpa3)
            
%             Fpa1b(logical([1,0,1; 0,1,0])) = Fpa2.get_data(NaN);
%             testCase.verifyEqual(Fpa1b, Fpa3)
            
            % =======================
            % Logical linear indexing
            % =======================
            % Index = 1D vector of logical.
            Fpa1a = Fpa([1,2,3; 4,5,6], NaN);
            Fpa1b = Fpa1a;
            Fpa2 = Fpa([9,8,7]', NaN);
            Fpa3 = Fpa([9,2,7; 4,8,6], NaN);
            Fpa1a(logical([1,0,0, 1,1,0]')) = Fpa2;
            testCase.verifyEqual(Fpa1a, Fpa3)

%             Fpa1b(logical([1,0,0, 1,1,0]')) = Fpa2.get_data(NaN);
%             testCase.verifyEqual(Fpa1b, Fpa3)
        end
        
        
        
        function test_subsasgn_different_size(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            % Copies FPAs, which requires FPA to be a non-handle class.
            assert(~isa(bicas.utils.FillPositionsArray.floatNan2logical([]), 'handle'))

            % Assign scalar to matrix.
            Fpa1a = Fpa([1,NaN,3; 4,5,NaN], NaN);
            Fpa1b = Fpa1a;
            Fpa2 = Fpa([9], NaN);
            Fpa3 = Fpa([9,9,9; 9,9,9], NaN);
            Fpa1a(:, :) = Fpa2;
            testCase.verifyEqual(Fpa1a, Fpa3)
%             Fpa1b(:, :) = Fpa2.get_data(NaN);
%             testCase.verifyEqual(Fpa1b, Fpa3)

        end


        function test_size_ndims(testCase)
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
                Fpa = bicas.utils.FillPositionsArray(v, 'FILL_VALUE', NaN);
                
                testCase.verifyEqual(isscalar(Fpa), isscalar(v))
                testCase.verifyEqual(size(Fpa),     size(v)    )
                testCase.verifyEqual(ndims(Fpa),    ndims(v)   )
                testCase.verifyEqual(isempty(Fpa),  isempty(v) )
                
                for iDim = 1:3
                    testCase.verifyEqual(size(Fpa, iDim), size(v, iDim) )
                end
            end    % for
        end
        
        
        
        function test_lt_gt_le_ge(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa

            Fpa1 = Fpa([1,  3, 2,   2, NaN], NaN);
            Fpa2 = Fpa([3,  1, 2, NaN,   2], NaN);
            FP = [false, false, false,  true,  true];

            ActFpa3 = Fpa1 < Fpa2;
            ExpFpa3 = bicas.utils.FillPositionsArray(...
                [true,  false, false, false, false], 'FILL_POSITIONS', FP);
            testCase.assertEqual(ActFpa3, ExpFpa3)

            ActFpa3 = Fpa1 > Fpa2;
            ExpFpa3 = bicas.utils.FillPositionsArray(...
                [false,  true, false, false, false], 'FILL_POSITIONS', FP);
            testCase.assertEqual(ActFpa3, ExpFpa3)

            ActFpa3 = Fpa1 <= Fpa2;
            ExpFpa3 = bicas.utils.FillPositionsArray(...
                [true,  false, true, false, false], 'FILL_POSITIONS', FP);
            testCase.assertEqual(ActFpa3, ExpFpa3)

            ActFpa3 = Fpa1 >= Fpa2;
            ExpFpa3 = bicas.utils.FillPositionsArray(...
                [false,  true, true, false, false], 'FILL_POSITIONS', FP);
            testCase.assertEqual(ActFpa3, ExpFpa3)
        end
        
        
        
        function test_floatNan2logical(testCase)
            function test_element_illegal_fail(mc2)
                Fpa = bicas.utils.FillPositionsArray.floatNan2logical(...
                    cast([2], mc2));
            end
            
            for mcCa = {'single', 'double'}'
                mc = mcCa{1};

                ExpFpa = bicas.utils.FillPositionsArray(...
                    [false, true,  false], 'FILL_POSITIONS', ...
                    [false, false,  true]);
                ActFpa = bicas.utils.FillPositionsArray.floatNan2logical(...
                    cast([0, 1, NaN], mc));
                testCase.verifyTrue(ExpFpa == ActFpa)

                testCase.verifyError(@() test_element_illegal_fail(mc), ?MException)
            end
        end
        
        
        
        function test_floatNan2int(testCase)
            for floatMcCa = {'single', 'double'}'                
                floatMc = floatMcCa{1};
                for intMcCa = {'uint8', 'int8'}'
                    intMc = intMcCa{1};

                    ExpFpa = bicas.utils.FillPositionsArray(...
                        cast([    0,     1,     0], intMc), 'FILL_POSITIONS', ...
                             [false, false,  true]);

                    ActFpa = bicas.utils.FillPositionsArray.floatNan2int(...
                        cast([0, 1, NaN], floatMc), intMc);
                    
                    testCase.verifyTrue(ExpFpa == ActFpa)
                end
            end
        end
        
        
        
        function test_cat_vertcat_horzcat(testCase)
            import bicas.utils.FillPositionsArray___UTEST.Fpa
            
            Fpa1 = Fpa([], nan);
            ActFpa = [Fpa1];
            ExpFpa = Fpa1;
            testCase.assertEqual(ActFpa, ExpFpa)
            
            Fpa1 = Fpa([1, 2], nan);
            Fpa2 = Fpa([3, 4], nan);

            ActFpa = [Fpa1, Fpa2];
            ExpFpa = Fpa([1, 2, 3, 4], nan);
            testCase.assertEqual(ActFpa, ExpFpa)

            ActFpa = [Fpa1; Fpa2];
            ExpFpa = Fpa([1, 2; 3, 4], nan);
            testCase.assertEqual(ActFpa, ExpFpa)

            ActFpa = cat(1, Fpa1, Fpa2);
            ExpFpa = Fpa([1, 2; 3, 4], nan);
            testCase.assertEqual(ActFpa, ExpFpa)

            ActFpa = cat(2, Fpa1, Fpa2);
            ExpFpa = Fpa([1, 2, 3, 4], nan);
            testCase.assertEqual(ActFpa, ExpFpa)
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        function Fpa = Fpa(dataAr, fv)
            Fpa = bicas.utils.FillPositionsArray(dataAr, 'FILL_VALUE', fv);
        end
        
        
        
        % Do multiple tests for FPAs which should be EQUAL.
        function test_equality(testCase, fpa1, fpa2, fv)
            % NOTE: For tests to be truly meaningful, implementation-dependent
            % fill position values should different.
            
            % Ensure that method "eq" is called.
            r = (fpa1 == fpa2);
            testCase.assertTrue(isscalar(r))
            testCase.verifyTrue(r)

            % Ensure that method "ne" is called.
            r = (fpa1 ~= fpa2);
            testCase.assertTrue(isscalar(r))
            testCase.verifyFalse(r)

            % Ensure that method "isequaln" is called.
            r = isequaln(fpa1, fpa2);
            testCase.assertTrue(isscalar(r))
            testCase.verifyTrue(r)

            % Verify equality through the methods for obtaining data.
            testCase.verifyEqual(...
                fpa1.get_data(fv), ...
                fpa2.get_data(fv))
            testCase.verifyEqual(...
                fpa1.fpAr, ...
                fpa2.fpAr)
        end
        


    end    % methods(Static, Access=private)



end
