%
% matlab.unittest automatic test code for
% solo.adm.get_directory_DSMDs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_directory_DSMDs___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_zero_direcctories(testCase)
            [ActDsmdArray, ActOiArray] = solo.adm.get_directory_DSMDs({});

            ExpOiArray = dir('~');
            ExpOiArray = ExpOiArray([], 1);

            testCase.assertEqual(ActDsmdArray, solo.adm.DSMD.empty(0, 1))
            testCase.assertEqual(ActOiArray, ExpOiArray)
        end



        function test_subdirectories(testCase)
            function subdirPath = create_subdir(parentDir, varargin)
                % PROPOSAL: Make into generic function.

                for i = 1:numel(varargin)
                    mkdir(parentDir, varargin{i})
                    parentDir = fullfile(parentDir, varargin{i});
                end

                subdirPath = parentDir;
            end

            function test(dirPathsCa, expDatasetPathsCa)
                [ActDsmdArray, ActOiArray] = solo.adm.get_directory_DSMDs(dirPathsCa);
                actPathCa = {ActDsmdArray.path}';

                testCase.assertEqual(...
                    sort(actPathCa(:)), ...
                    sort(expDatasetPathsCa(:)))
                testCase.assertEqual(...
                    size(ActOiArray), size(ActDsmdArray))
            end

            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture)
            testDir = pwd;
            cd('~')    % Move to any OTHER unrelated directory.

            junk = create_subdir(testDir, 'subdir1', 'subdir11');
            dir2 = create_subdir(testDir, 'subdir2', 'subdir21');
            dir3 = create_subdir(testDir, 'subdir3', 'subdir31');

            % NOTE: No file in subdir1.
            file2  = fullfile(dir2, 'solo_L1R_rpw-lfr-surv-cwf-e_20240101_V02.cdf');
            file3a = fullfile(dir3, 'solo_L1R_rpw-lfr-surv-cwf-e_20240102_V03.cdf');
            file3b = fullfile(dir3, 'solo_L1R_rpw-lfr-surv-cwf-e_20240103_V04.cdf');
            file3c = fullfile(dir3, 'NOT_DATASET.DAT');
            irf.fs.create_empty_file(file2)
            irf.fs.create_empty_file(file3a)
            irf.fs.create_empty_file(file3b)
            irf.fs.create_empty_file(file3c)

            test(...
                {testDir}, ...
                {file2, file3a, file3b})

            test(...
                {fullfile(testDir, 'subdir1')}, ...
                {})
            test(...
                {fullfile(testDir, 'subdir2')}, ...
                {file2})

            test(...
                {fullfile(testDir, 'subdir3')}, ...
                {file3a, file3b})
        end



    end    % methods(Test)



end
