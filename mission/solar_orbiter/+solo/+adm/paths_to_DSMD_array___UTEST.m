%
% matlab.unittest automatic test code for
% solo.adm.paths_to_DSMD_array().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef paths_to_DSMD_array___UTEST < matlab.unittest.TestCase

    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)
            % NOTE: Must test both with and without
            %   (1) "-cdag" and/or (2) parent directories in the path.
            %
            % Example datasets from flight:
            %   solo_HK_rpw-bia_20200301_V01.cdf                   # NOTE: No -cdag.
            %   solo_L2_rpw-lfr-surv-cwf-e-cdag_20200213_V01.cdf   # NOTE: -cdag.
            %   solo_L1_rpw-bia-sweep-cdag_20200307T053018-20200307T053330_V01.cdf
            %   solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf

            function test(filePathList, expDsmdArray)
                actDsmdArray = solo.adm.paths_to_DSMD_array(filePathList);
                testCase.assertEqual(expDsmdArray, actDsmdArray)
            end
            
            function dt = dtu(varargin)
                 dt = datetime(varargin{:}, 'TimeZone', 'UTCLeapSeconds');
            end

            DSMD_1 = solo.adm.DSMD(...
                'solo_HK_rpw-bia_20200301_V01.cdf', ...
                'SOLO_HK_RPW-BIA', 1, false, ...
                dtu(2020, 03, 01, 00, 00, 00), ...
                dtu(2020, 03, 02, 00, 00, 00));

            DSMD_2 = solo.adm.DSMD(...
                'dir/solo_L1_rpw-bia-sweep-cdag_20200307T053018-20200307T053330_V02.cdf', ...
                'SOLO_L1_RPW-BIA-SWEEP', 2, true, ...
                dtu(2020, 03, 07, 05, 30, 18), ...
                dtu(2020, 03, 07, 05, 33, 30));

            % New filenaming convention for currents. EJ+XB-email 2020-05-27.
            DSMD_3 = solo.adm.DSMD(...
                'solo_L1_rpw-bia-current-cdag_20200301-20200331_V03.cdf', ...
                'SOLO_L1_RPW-BIA-CURRENT', 3, true, ...
                dtu(2020, 03, 01, 00, 00, 00), ...
                dtu(2020, 04, 01, 00, 00, 00));

            FILE_CDF_IGNORE    = '/temp/BIAS_RCT.cdf';
            FILE_NONCDF_IGNORE = '/temp/readme.txt';



            % Empty argument (no files)
            test({}, solo.adm.DSMD.empty(0,1));

            test({FILE_CDF_IGNORE},    solo.adm.DSMD.empty(0,1));
            test({FILE_NONCDF_IGNORE}, solo.adm.DSMD.empty(0,1));

            test(...
                {DSMD_1.path}, ...
                [DSMD_1]);

            test(...
                {DSMD_2.path}, ...
                [DSMD_2]);

            test(...
                {DSMD_3.path}, ...
                [DSMD_3]);

            test(...
                {DSMD_1.path; FILE_NONCDF_IGNORE; DSMD_2.path; FILE_CDF_IGNORE; DSMD_3.path}, ...
                [DSMD_1; DSMD_2; DSMD_3]);
        end



    end    % methods(Test)

end
