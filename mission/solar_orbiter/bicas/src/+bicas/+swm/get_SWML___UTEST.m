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
                Bso = bicas.swm.get_SWML___UTEST.get_BSO(...
                    l1l2Enabled, l2l2Enabled, l2l3Enabled);

                Swml = bicas.swm.get_SWML(Bso);

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

        function Bso = get_BSO(l1l2Enabled, l2l2Enabled, l2l3Enabled)
            Bso = bicas.create_default_BSO();

            Bso.override_value('SWM.L1-L2_ENABLED',         l1l2Enabled, 'test')
            Bso.override_value('SWM.L2-L2_CWF-DSR_ENABLED', l2l2Enabled, 'test')
            Bso.override_value('SWM.L2-L3_ENABLED',         l2l3Enabled, 'test')

            Bso.make_read_only()
        end

    end    % methods(Static, Access=private)



end
