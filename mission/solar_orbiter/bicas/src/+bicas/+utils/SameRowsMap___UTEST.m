%
% matlab.unittest automatic test code for bicas.utils.SameRowsMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SameRowsMap___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test_basic(testCase)
            
            % Adds key-values, zero rows, char keys, 'empty'.
            V1 = zeros(0, 1);
            V2 = ones( 0, 1, 2);
            
            M = bicas.utils.SameRowsMap('char', 0, 'empty');
            testCase.assertEqual(M.keys(), cell(0, 1))
            testCase.assertEqual(M.length, 0)
            testCase.assertEqual(M.nRows,  0)
            
            M.add('K1', V1)
            testCase.assertEqual(M.keys,   {'K1'})
            testCase.assertEqual(M.length, 1)
            
            testCase.assertFalse(M.isKey('K2'))
            M.add('K2', V2)
            testCase.assertTrue(bicas.utils.SameRowsMap.key_sets_equal(...
                M.keys(), {'K1', 'K2'} ...
            ))
            testCase.assertEqual(M.length, 2)
            testCase.assertTrue(M.isKey('K2'))
            
            testCase.assertEqual(M.get('K1'), V1);
            testCase.assertEqual(M.get('K2'), V2);
            testCase.assertEqual(M.nRows,  0)
            

            
            % Zero number of constant values.
            M = bicas.utils.SameRowsMap('double', 3, 'constant', [1;2;3], {});
            testCase.assertEqual(M.nRows, 3)
            

            
            % double keys, non-zero rows, 'constant'
            M = bicas.utils.SameRowsMap('double', 3, 'constant', [1;2;3], {9});
            testCase.assertEqual(M.nRows, 3)
            testCase.assertEqual(M.get(9), [1;2;3])
            
            

            % Initial value has inconsistent number of rows.
            testCase.assertError(...
                @() (bicas.utils.SameRowsMap('double', 3, 'constant', [1;2])), ...
                ?MException)
        end
        
        
        
        function test_setRows(testCase)
            
            % Insert zero rows into zero rows.
            M1 = bicas.utils.SameRowsMap('char', 0, 'empty');
            M1.add('K2', zeros(0,0))
            M2 = bicas.utils.SameRowsMap('char', 0, 'empty');
            M2.add('K2', zeros(0,0))
            M1.setRows(M2, zeros(0,1))
            testCase.assertEqual(M1.get('K2'), zeros(0,0))

            % Insert zero rows into non-zero rows.
            M1 = bicas.utils.SameRowsMap('char', 3, 'empty');
            M1.add('K2', zeros(3,0))
            M2 = bicas.utils.SameRowsMap('char', 0, 'empty');
            M2.add('K2', zeros(0,0))
            M1.setRows(M2, zeros(0,1))
            testCase.assertEqual(M1.get('K2'), zeros(3,0))

            % Preserve type
            M1 = bicas.utils.SameRowsMap('char', 4, 'empty');
            M1.add('K2', int16([1;2;3;4]))
            M2 = bicas.utils.SameRowsMap('char', 2, 'empty');
            M2.add('K2', int16([[8;9]]))
            M1.setRows(M2, [2;3])
            testCase.assertEqual(M1.get('K2'), int16([1;8;9;4]))



            % =======================
            % Higher-dimensional case
            % =======================
            % IMPLEMENTATION NOTE: Simpler to just define one example variable
            % (of higher dimensionality) and then extract smaller variables from
            % it.
            V1(:, :, 1) = [...
                111, 121, 131; ...
                211, 221, 231; ...
                311, 321, 331; ...
                411, 421, 431 ...
            ];
            V1(:, :, 2) = [...
                112, 122, 132; ...
                212, 222, 232; ...
                312, 322, 332; ...
                412, 422, 432 ...
            ];
        
            M1 = bicas.utils.SameRowsMap('char', 4, 'constant', V1,            {'K'});
            M2 = bicas.utils.SameRowsMap('char', 2, 'constant', V1(1:2, :, :), {'K'});
            M1.setRows(M2, [3;2])    % NOTE: Decrementing indices.
            
            V2              = V1;
            V2([3,2], :, :) = V1(1:2, :, :);
            testCase.assertEqual(M1.get('K'), V2)



            % =======
            % Illegal
            % =======
            % Incompatible types
            M1 = bicas.utils.SameRowsMap('char', 2, 'constant', int8([1;2]), {'K'});
            M2 = bicas.utils.SameRowsMap('char', 1, 'constant', [9],         {'K'});
            testCase.assertError(...
                @() (M1.setRows(M2, [1])), ...
                ?MException)
            % Incompatible array sizes
            M1 = bicas.utils.SameRowsMap('char', 2, 'constant', [1 2;3 4], {'K'});
            M2 = bicas.utils.SameRowsMap('char', 1, 'constant', [9],       {'K'});
            testCase.assertError(...
                @() (M1.setRows(M2, [1])), ...
                ?MException)
            % Different sets of keys.
            M1 = bicas.utils.SameRowsMap('char', 2, 'constant', [1;2], {'K1', 'K2a'});
            M2 = bicas.utils.SameRowsMap('char', 1, 'constant', [9],   {'K1', 'K2b'});
            testCase.assertError(...
                @() (M1.setRows(M2, [1])), ...
                ?MException)
        end
        
        
        
        function test_key_sets_equal(testCase)
            
            % One output variable.
            function test(keysCa1, keysCa2, expEqual)
                actEqual = bicas.utils.SameRowsMap.key_sets_equal(keysCa1, keysCa2);
                testCase.assertEqual(actEqual, expEqual)
            end

            test({}, {}, true)
            test({ 1 }, { 1 }, true)
            test({'1'}, {'1'}, true)
            test({ 1 }, { 2 }, false)
            test({'1'}, {'2'}, false)
            
            test({'asd', 1}, {'asd', 1}, true)
            test({'asd', 1}, {'asd', 2}, false)
            test({'asd', 1}, {'ASD', 1}, false)
        end
        
        
        
        % NOTE: Does not test method per se.
        function test_equal(testCase)
            % TODO: NaN
            
            M1a = bicas.utils.SameRowsMap('char',   1, 'empty');
            M1b = bicas.utils.SameRowsMap('char',   1, 'empty');
            M2  = bicas.utils.SameRowsMap('char',   2, 'empty');
            M3  = bicas.utils.SameRowsMap('double', 1, 'empty');
            
            testCase.assertTrue( M1a == M1b)
            testCase.assertFalse(M1a == M2)
            testCase.assertFalse(M1a == M3)
            
            M1a = bicas.utils.SameRowsMap('char', 1, 'constant', [1],   {'K1', 'K2'});
            M1b = bicas.utils.SameRowsMap('char', 1, 'constant', [1],   {'K1', 'K2'});
            M2  = bicas.utils.SameRowsMap('char', 1, 'constant', [9],   {'K1', 'K2'});
            M3  = bicas.utils.SameRowsMap('char', 2, 'constant', [1;2], {'K1', 'K2'});
            
            testCase.assertTrue( (M1a == M1b))
            testCase.assertFalse((M1a == M2))
            testCase.assertFalse((M1a == M3))
            
            % NaN, different key order.
            M1a = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1a.add('K1', [NaN])
            M1a.add('K2', [2])
            M1b = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1b.add('K2', [2])
            M1b.add('K1', [NaN])
            
            testCase.assertTrue( M1a == M1b)
            testCase.assertTrue( M1a == M1b)
            testCase.assertFalse(M1a == M2 )
            
            % Different MATLAB classes.
            M1  = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1.add('K1', [1])
            M2  = bicas.utils.SameRowsMap('char', 1, 'empty');
            M2.add('K1', int8([1]))
            testCase.assertFalse(M1a == M2 )
        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
    end    % methods(Static, Access=private)

    
    
end
