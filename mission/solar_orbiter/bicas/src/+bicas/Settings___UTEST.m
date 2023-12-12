%
% matlab.unittest automatic test code for bicas.Settings class.
%
% NOTE: Does not test cell-valued settings since they seem not to be fully
%       supported. To be investigated. /2023-12-12
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Settings___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_empty(testCase)
            S = bicas.Settings();
            S.disable_define()
            S.make_read_only()
        end



        % Test that .get_fv() can not access values when it is not supposed
        % to (should raise error).
        function test_get_fv_access(testCase)
            KEY   = 'TEST_KEY';
            VALUE = 99;

            S = bicas.Settings();
            S.define_setting(KEY, VALUE);

            testCase.verifyError(...
                @() S.get_fv(KEY), ...
                ?MException)

            S.disable_define()

            testCase.verifyError(...
                @() S.get_fv(KEY), ...
                ?MException)

            S.make_read_only()

            actValue = S.get_fv(KEY);
            testCase.assertEqual(actValue, VALUE);
        end



        % Test storing different value types, with and without overriding.
        %
        % Tests multiple methods:
        %   .override_value()
        %   .override_values_from_strings()
        %   .get_fv()
        %   .get_final_value_array()
        %   etc.
        %
        % NOTE: Does not test cell arrays.
        function test_value_types(testCase)

            % (1) set a key-value,
            % (not override value)
            % (2) read value.
            function test_initial(initialValue)
                KEY = 'TEST_KEY';

                S = bicas.Settings();

                S.define_setting(KEY, initialValue)
                S.disable_define()
                S.make_read_only()

                actValue = S.get_fv(KEY);
                testCase.assertEqual(actValue, initialValue)
            end

            % (1) set a key-value,
            % (2) override value,
            % (3) read value.
            %
            % ARGUMENTS
            % =========
            % overrideValue, overrideValueAsString
            %       Non-string representation and string representation of the
            %       same value.
            function test_override(initialValue, overrideValue, overrideValueAsString)
                assert(ischar(overrideValueAsString))

                KEY = 'TEST_KEY';

                function S = pre_override()
                    S = bicas.Settings();

                    S.define_setting(KEY, initialValue)
                    S.disable_define()
                end

                function post_override(S)
                    actValue = S.get_fv(KEY);
                    testCase.assertEqual(actValue, overrideValue)

                    actS = S.get_final_value_array(KEY);
                    expS = struct(...
                        'value',       {initialValue, overrideValue}, ...
                        'valueSource', {'default',    'test'});
                    testCase.assertEqual(actS, expS)
                end

                %==================================
                % Override using .override_value()
                %==================================
                if 1
                    S1 = pre_override();

                    S1.override_value(KEY, overrideValue, 'test')
                    S1.make_read_only()

                    post_override(S1);
                end

                %================================================
                % Override using .override_values_from_strings()
                %================================================
                if 1
                    S2 = pre_override();

                    ModifiedSettingsMap      = containers.Map();
                    ModifiedSettingsMap(KEY) = overrideValueAsString;
                    S2.override_values_from_strings(ModifiedSettingsMap, 'test');
                    S2.make_read_only()

                    % NOTE: Checks that overriding using "overrideValueAsString"
                    % yields "overrideValue", i.e. use both arguments which must
                    % thus represent the same value.
                    post_override(S2);
                end
            end

            %====================
            % Test default value
            %====================
            if 1
                test_initial(false)
                test_initial(true)
                test_initial('')
                test_initial('Abc')
                % NOTE: 0x0 array is not allowed.
                test_initial(ones(1,0))
                test_initial(ones(0,1))
                test_initial(123)
                test_initial([123, 234, 345])
            end

            %=====================
            % Test override value
            %=====================
            test_override(true,      false,      'false')
            test_override(false,     true,       'true')
            test_override('',        'Abc',      'Abc')
            test_override('initial', '',         '')
            test_override('initial', 'Abc',      'Abc')

            test_override(ones(1,0), [2], '2')    % Empty initial array.

            % NOTE: Empty string (0x0 char) interpreted as numeric value is
            % converted to empty row vector (of numeric).
            test_override(99,        ones(1,0),  '')
            test_override(99,        123,        '123')
            test_override(99,        [123, 234], '123, 234')
            test_override(99,        [123, 234], '[123, 234]')
        end



        function test_get_keys(testCase)
            %===========
            % Zero keys
            %===========
            S1 = bicas.Settings();
            testCase.assertEqual(S1.get_keys(), cell(1,0))

            S1.disable_define()
            testCase.assertEqual(S1.get_keys(), cell(1,0))

            S1.make_read_only()
            testCase.assertEqual(S1.get_keys(), cell(1,0))

            %===============
            % Non-zero keys
            %===============
            KEYS_SET = {'KEY_1', 'KEY_2'};
            S2 = bicas.Settings();
            S2.define_setting(KEYS_SET{1}, 99)
            S2.define_setting(KEYS_SET{2}, 'Abc')
            testCase.assertEqual(S2.get_keys(), KEYS_SET)

            S2.disable_define()
            testCase.assertEqual(S2.get_keys(), KEYS_SET)

            S2.make_read_only()
            testCase.assertEqual(S2.get_keys(), KEYS_SET)
        end



        % Test .get_final_value_array() for non-overridden values.
        %
        % NOTE: test_value_types() also tests the same method but for overridden
        % values.
        function test_get_final_value_array(testCase)
            KEY_1 = 'KEY_1';
            KEY_2 = 'KEY_2';
            VALUE_1 = 99;
            VALUE_2 = 'Abc';

            S1 = bicas.Settings();
            S1.define_setting(KEY_1, VALUE_1);
            S1.make_read_only()

            actS = S1.get_final_value_array(KEY_1);
            expS = struct('value', VALUE_1, 'valueSource', 'default');
            testCase.assertEqual(actS, expS)

            % Two simultaneous keys

            S2 = bicas.Settings();
            S2.define_setting(KEY_1, VALUE_1);
            S2.define_setting(KEY_2, VALUE_2);
            S2.make_read_only()

            actS = S2.get_final_value_array(KEY_2);
            expS = struct('value', VALUE_2, 'valueSource', 'default');
            testCase.assertEqual(actS, expS)

            actS = S2.get_final_value_array(KEY_1);
            expS = struct('value', VALUE_1, 'valueSource', 'default');
            testCase.assertEqual(actS, expS)
        end



    end    % methods(Test)



end
