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



        function test_misc(testCase)
            import bicas.utils.FillPositionsArray___UTEST.test_equality

            % ===========================================================
            % Test legal & illegal read & write to *properties* (public &
            % private)
            % ===========================================================
            Fpa = bicas.utils.FillPositionsArray([], 'fill value', 0);
            
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
                % Does not work, since class does not work.
                testCase.verifyError(@() (Fpa.dataAr), ?MException)
            end
            
            
            %=============
            % 0x0, double
            %=============
            if 1
                Fpa1 = bicas.utils.FillPositionsArray([], 'fill value', 0);
                Fpa2 = bicas.utils.FillPositionsArray([], 'fill positions', logical([]));
                %
                test_equality(testCase, Fpa1, Fpa2, -99)
                testCase.verifyEqual(Fpa1.get_data(999), [])
                testCase.verifyEqual(Fpa1.fpAr,          logical([]))
            end

            %===========
            % 1D, int64
            %===========
            if 1
                Fpa1 = bicas.utils.FillPositionsArray(...
                    int64([1, 2, 3, -9]), 'fill value', int64(-9));
                Fpa2 = bicas.utils.FillPositionsArray(...
                    int64([1, 2, 3, -9]), 'fill positions', logical([0, 0, 0, 1]));
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
            if 1
                Fpa1 = bicas.utils.FillPositionsArray(...
                    ['abX'; 'dXf'], 'fill value', 'X');
                Fpa2 = bicas.utils.FillPositionsArray(...
                    ['abY'; 'dYf'], 'fill positions', logical([0, 0, 1; 0, 1, 0]));
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
                    [1, 2, fv1; 4, fv1, 6], 'fill value', fv1);
                Fpa2 = bicas.utils.FillPositionsArray(...
                    [1, 2,  -2; 4,  -5, 6], 'fill positions', ...
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
                                    [x, fv1], 'fill value', fv1);
                                Fpa2 = bicas.utils.FillPositionsArray(...
                                    [x,   x], 'fill positions', ...
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



        function test_convert(testCase)
            % Convert to type that forbids negative values, while input used
            % negative FV.
            Fpa1 = bicas.utils.FillPositionsArray([0,1,-1,2,3], 'fill value', -1);
            Fpa2   = Fpa1.convert(@(x) (x), 'uint16', 99);
            dataAr = Fpa2.get_data(uint16(999));
            testCase.verifyEqual(dataAr, uint16([0,1,999,2,3]))
            
            % Operation that raises error for NaN.
            % Convert to type that forbids NaN, while input used NaN as FV.
            Fpa1 = bicas.utils.FillPositionsArray([0, 3, NaN], 'fill value', NaN);
            Fpa2   = Fpa1.convert(@(x) (~x), 'logical', 99);
            dataAr = Fpa2.get_data(true);
            testCase.verifyEqual(dataAr, [true, false, true])
        end
        
        
        
        function test_cast(testCase)
            % float FPA --> logical FPA
            Fpa1   = bicas.utils.FillPositionsArray([0, 1, NaN], 'fill value', NaN);
            Fpa2   = Fpa1.cast('logical', 0);
            dataAr = Fpa2.get_data(true);
            testCase.verifyEqual(dataAr, [false, true, true])
            
            % logical FPA --> float FPA
            Fpa1   = bicas.utils.FillPositionsArray([false, true, false], 'fill positions', [false, false, true]);
            Fpa2   = Fpa1.cast('double', false);
            dataAr = Fpa2.get_data(NaN);
            testCase.verifyEqual(dataAr, [0, 1, NaN])
        end
        
        
        
        function test_subsref(testCase)
            % ===
            % 0x0
            % ===
            Fpa1 = bicas.utils.FillPositionsArray([], 'fill value', NaN);
            Fpa2 = Fpa1(:);    % 0x0 --> 0x1
            dataAr2 = Fpa2.get_data(NaN);
            testCase.verifyEqual(dataAr2, ones(0,1))
            
            % ==
            % 2D
            % ==
            Fpa1 = bicas.utils.FillPositionsArray([1,2; 3,4; 5,NaN], 'fill value', NaN);
            
            Fpa2 = Fpa1(3, 1:2);    % 2 values, one range.
            dataAr2 = Fpa2.get_data(NaN);
            testCase.verifyEqual(dataAr2, [5, NaN])
            
            Fpa3 = Fpa1(:);    % 1D indexing.
            dataAr3 = Fpa3.get_data(NaN);
            testCase.verifyEqual(dataAr3, [1,3,5,2,4,NaN]')
        end
        
        
        
        function test_size(testCase)
            % IMPLEMENTATION NOTE: In a sense, this does not only test the code,
            % but also the author's understanding of "overloading" with a method
            % size().
            
            TEST_DATA_CA = {...
                ones(0, 0), ...
                ones(1, 1), ...
                ones(1, 1, 0), ...
                ones(3, 1), ...
                ones(1, 3), ...
                ones(4, 5), ...
            };
            
            for ca = TEST_DATA_CA(:)'
                v = ca{1};
                Fpa = bicas.utils.FillPositionsArray(v, 'fill value', NaN);
                
                testCase.verifyEqual(isscalar(Fpa), isscalar(v))
                testCase.verifyEqual(size(Fpa),     size(v)    )
                testCase.verifyEqual(ndims(Fpa),    ndims(v)   )
                
                for iDim = 1:3
                    testCase.verifyEqual(size(Fpa, iDim),  size(v, iDim) )
                end
            end
            
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Do multiple tests for FPAs which should be equal.
        function test_equality(testCase, fpa1, fpa2, fillValue)
            % Ensure that method "eq" is called.
            r = (fpa1 == fpa2);
            testCase.assertTrue(isscalar(r))
            testCase.verifyTrue(r)

            % Ensure that method "ne" is called.
            r = (fpa1 ~= fpa2);
            testCase.assertTrue(isscalar(r))
            testCase.verifyFalse(r)

            testCase.verifyEqual(...
                fpa1.get_data(fillValue), ...
                fpa2.get_data(fillValue))
            testCase.verifyEqual(...
                fpa1.fpAr, ...
                fpa2.fpAr)
        end
        


    end    % methods(Static, Access=private)



end
