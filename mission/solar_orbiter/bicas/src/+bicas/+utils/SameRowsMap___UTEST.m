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
            % Test sequences of operations on a single Map.
            % Exlude setRows().
            
            % ==============================================
            % Adds key-values, zero rows, char keys, 'empty'
            % ==============================================
            V1 = zeros(0, 1);
            V2 = ones( 0, 1, 2);
            
            M = bicas.utils.SameRowsMap('char', 0, 'empty');
            
            testCase.assertEqual(M.keys(), cell(0, 1))
            testCase.assertEqual(M.length, 0)
            testCase.assertEqual(M.nRows,  0)
            bicas.utils.SameRowsMap___UTEST.test_keys_values(testCase, M, {}, {})

            M.add('K1', V1)
            
            testCase.assertEqual(M.keys,   {'K1'})
            testCase.assertEqual(M.length, 1)
            bicas.utils.SameRowsMap___UTEST.test_keys_values(testCase, M, {'K1'}, {V1})            
            testCase.assertFalse(M.isKey('K2'))
            
            M.add('K2', V2)
            
            testCase.assertEqual(M.length, 2)
            testCase.assertTrue( M.isKey('K2'))
            
            % Test different orders. ==> Effectively testing the helper
            % function "test_keys_values()".
            bicas.utils.SameRowsMap___UTEST.test_keys_values(testCase, M, {'K2', 'K1'}, {V2, V1})
            bicas.utils.SameRowsMap___UTEST.test_keys_values(testCase, M, {'K1', 'K2'}, {V1, V2})
            
            testCase.assertEqual(M('K1'), V1);
            testCase.assertEqual(M('K2'), V2);
            testCase.assertEqual(M.nRows,  0)
            
            % =================================================
            % Zero number of constant values (test constructor)
            % =================================================
            M = bicas.utils.SameRowsMap('double', 3, 'constant', [1;2;3], {});
            testCase.assertEqual(M.nRows, 3)

            % ======================================
            % double keys, non-zero rows, 'constant'
            % ======================================
            V = [1;2;3];
            M = bicas.utils.SameRowsMap('double', 3, 'constant', V, {9});
            testCase.assertEqual(M.nRows, 3)
            testCase.assertEqual(M(9), V)
            bicas.utils.SameRowsMap___UTEST.test_keys_values(testCase, M, {9}, {V})

            % =============================================
            % Initial value has inconsistent number of rows
            % =============================================
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
            testCase.assertEqual(M1('K2'), zeros(0,0))
            
            % Insert zero rows into non-zero rows.
            M1 = bicas.utils.SameRowsMap('char', 3, 'empty');
            M1.add('K2', zeros(3,0))
            M2 = bicas.utils.SameRowsMap('char', 0, 'empty');
            M2.add('K2', zeros(0,0))
            M1.setRows(M2, zeros(0,1))
            testCase.assertEqual(M1('K2'), zeros(3,0))

            % Preserve type
            M1 = bicas.utils.SameRowsMap('char', 4, 'empty');
            M1.add('K2', int16([1;2;3;4]))
            M2 = bicas.utils.SameRowsMap('char', 2, 'empty');
            M2.add('K2', int16([[8;9]]))
            M1.setRows(M2, [2;3])
            testCase.assertEqual(M1('K2'), int16([1;8;9;4]))



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
            testCase.assertEqual(M1('K'), V2)



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
        
        
        
        % NOTE: Does not test method per se.
        function test_equal(testCase)
            
            % NOTE: Having no keys with different KeyType means that the key
            % values do not distinguish the KeyTypes.
            M1a = bicas.utils.SameRowsMap('char',   1, 'empty');
            M1b = bicas.utils.SameRowsMap('char',   1, 'empty');
            M2a = bicas.utils.SameRowsMap('double', 1, 'empty');
            M2b = bicas.utils.SameRowsMap('double', 1, 'empty');
            M3  = bicas.utils.SameRowsMap('char',   2, 'empty');
            
            testCase.assertTrue( M1a == M1b)
            testCase.assertFalse(M1a == M2a)
            testCase.assertTrue( M2a == M2b)
            testCase.assertFalse(M1a == M3 )
            
            M1a = bicas.utils.SameRowsMap('char',   1, 'constant', [1],   {'K1', 'K2'});
            M1b = bicas.utils.SameRowsMap('char',   1, 'constant', [1],   {'K1', 'K2'});
            M2a = bicas.utils.SameRowsMap('double', 1, 'constant', [9],   {1, 2});
            M2b = bicas.utils.SameRowsMap('double', 1, 'constant', [9],   {1, 2});            
            M3  = bicas.utils.SameRowsMap('char',   1, 'constant', [9],   {'K1', 'K2'});
            M4  = bicas.utils.SameRowsMap('char',   2, 'constant', [1;2], {'K1', 'K2'});

            testCase.assertTrue( M1a == M1b)
            testCase.assertTrue( M2a == M2b)
            testCase.assertFalse(M1a == M3 )
            testCase.assertFalse(M1a == M4)
            testCase.assertFalse(M2a == M3)
            testCase.assertFalse(M2a == M4)
            
            % NaN, different key order.
            M1a = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1a.add('K1', [NaN])
            M1a.add('K2', [2])
            M1b = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1b.add('K2', [2])
            M1b.add('K1', [NaN])
            
            testCase.assertTrue( M1a == M1b)
            testCase.assertTrue( M1a == M1b)
            testCase.assertFalse(M1a == M3 )
            
            % Different MATLAB classes.
            M1  = bicas.utils.SameRowsMap('char', 1, 'empty');
            M1.add('K1', [1])
            M3  = bicas.utils.SameRowsMap('char', 1, 'empty');
            M3.add('K1', int8([1]))
            testCase.assertFalse(M1a == M3 )
        end
        
        
        
        function test_display(testCase)
            M  = bicas.utils.SameRowsMap('char', 0, 'constant', uint8(ones(0,2)), {'K1'});
            disp(M)
            
            M  = bicas.utils.SameRowsMap('double', 1, 'constant', [9], {1});
            disp(M)
            
            M  = bicas.utils.SameRowsMap('char', 3, 'constant', [1,2;3,4;5,6], {'K1', 'K2'});
            disp(M)
        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Helper function. Test that Map returns correct values for methods
        % .keys() and .values(). Takes into account that keys/values may be in
        % a different order than specified.
        function test_keys_values(testCase, Map, expKeysCa, expValuesCa)
            % NOTE: Implementation should permit both numbers and strings as
            % keys, mixed.
            
            actKeysCa   = Map.keys();
            actValuesCa = Map.values();
            
            testCase.assertTrue(bicas.utils.object_sets_isequaln(actKeysCa, expKeysCa))
            
            for iAct = 1:numel(actKeysCa)
                for iExp = 1:numel(expKeysCa)
                    if isequal(actKeysCa{iAct}, expKeysCa{iExp})
                        actValue = actValuesCa{iAct};
                        expValue = expValuesCa{iExp};
                        testCase.assertEqual(actValue, expValue)
                    end
                end
            end
            
        end
        
        
        
    end    % methods(Static, Access=private)

    
    
end
