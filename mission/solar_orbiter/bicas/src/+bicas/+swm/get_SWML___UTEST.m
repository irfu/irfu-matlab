%
% matlab.unittest automatic test code for bicas.swm.get_SWML().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_SWML___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            function test(l1l2Enabled, l2l2Enabled, l2l3Enabled)
                % Tests (1) non-crash, (2) class of return value.
                L = bicas.Logger('none', false);
                SETTINGS = bicas.swm.get_SWML___UTEST.get_SETTINGS(...
                    l1l2Enabled, l2l2Enabled, l2l3Enabled);

                Swml = bicas.swm.get_SWML(SETTINGS, L);

                testCase.verifyClass(Swml, 'bicas.swm.SoftwareModeList')
            end

            %===================================================================

            for l1l2Enabled = [0, 1]
                for l2l2Enabled = [0, 1]
                    for l2l3Enabled = [0, 1]
                        test(l1l2Enabled, l2l2Enabled, l2l3Enabled)
                    end
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

        function SETTINGS = get_SETTINGS(l1l2Enabled, l2l2Enabled, l2l3Enabled)
            SETTINGS = bicas.create_default_SETTINGS();

            SETTINGS.override_value('SWM.L1-L2_ENABLED',         l1l2Enabled, 'test')
            SETTINGS.override_value('SWM.L2-L2_CWF-DSR_ENABLED', l2l2Enabled, 'test')
            SETTINGS.override_value('SWM.L2-L3_ENABLED',         l2l3Enabled, 'test')

            SETTINGS.make_read_only()
        end

    end    % methods(Static, Access=private)



end
