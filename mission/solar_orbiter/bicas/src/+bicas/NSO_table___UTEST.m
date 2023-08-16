%
% matlab.unittest automatic test code for bicas.NSO_table.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef NSO_table___UTEST < matlab.unittest.TestCase
    % PROPOSAL: Extend with more tests.
    %   PROPOSAL: Empty (legal) NSO file?
    %       NOTE: Already test-loading default NSO XML file.



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        % % Test constructor.
        % % NOTE: test_get_NSO_timestamps() indirectly tests the constructor.
        % function test_NSO_table(testCase)
        % end



        % Test method bicas.NSO_table.get_NSO_timestamps(), but indirectly
        % uses/tests bicas.NSO_table() (constructor) by its nature.
        function test_get_NSO_timestamps(testCase)
            % PROBLEM: How handle that return value may change the order of
            %          events depending on implementation?

            function test(...
                    evtStartTt2000Array, evtStopTt2000Array, evtNsoIdCa, tt2000Array, ...
                    expBEvtArraysCa, expEvtNsoIdCa, expIGlobalEventsArray)

                % Normalize input
                evtStartTt2000Array = int64(evtStartTt2000Array(:));
                evtStopTt2000Array  = int64(evtStopTt2000Array(:));
                evtNsoIdCa          = evtNsoIdCa(:);
                tt2000Array         = int64(tt2000Array(:));
                % Normalize output
                expBEvtArraysCa       = expBEvtArraysCa(:);
                expEvtNsoIdCa         = expEvtNsoIdCa(:);
                expIGlobalEventsArray = expIGlobalEventsArray(:);

                NsoTable = bicas.NSO_table(evtStartTt2000Array, evtStopTt2000Array, evtNsoIdCa);

                [actBEvtArraysCa, actEvtNsoIdCa, actIGlobalEventsArray] = NsoTable.get_NSO_timestamps(tt2000Array);
                testCase.verifyEqual(actBEvtArraysCa,       expBEvtArraysCa)
                testCase.verifyEqual(actEvtNsoIdCa,         expEvtNsoIdCa)
                testCase.verifyEqual(actIGlobalEventsArray, expIGlobalEventsArray)
            end

            %===================================================================

            NSOID_1 = bicas.constants.NSOID.FULL_SATURATION;
            NSOID_2 = bicas.constants.NSOID.PARTIAL_SATURATION;

            ENA = zeros(0, 1);   % Empty Numeric Array
            ECA = cell(0, 1);    % Empty Cell    Array

            % Test every combination of
            % (1) empty & (2) non-empty
            % for (a) NSO table & (b) submitted timestamps
            % without any overlap (empty output).
            for tt2000ArrayCa = {ENA, [100:200]'}
                tt2000Array = int64(tt2000ArrayCa{1});

                % Empty NSO table.
                test(...
                    ENA, ENA, ECA, ...
                    tt2000Array, ...
                    ECA, ECA, ENA)

                % Non-empty NSO table.
                test(...
                    [10], [20], {NSOID_1}, ...
                    tt2000Array, ...
                    ECA, ECA, ENA)
            end

            % NSOs do not overlap beginning & end.
            test(...
                [1, 5], [2, 7], {NSOID_1, NSOID_2}, ...
                [0:9], ...
                {...
                    logical([0, 1, 1, 0, 0, 0, 0, 0, 0, 0]'), ...
                    logical([0, 0, 0, 0, 0, 1, 1, 1, 0, 0]') ...
                }, ...
                {NSOID_1, NSOID_2}, [1, 2])

            % NSOs overlap beginning & end.
            test(...
                [-1, 5], [2, 12], {NSOID_1, NSOID_2}, ...
                [0:9], ...
                {...
                    logical([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]'), ...
                    logical([0, 0, 0, 0, 0, 1, 1, 1, 1, 1]') ...
                }, ...
                {NSOID_1, NSOID_2}, [1, 2])

            % Overlapping NSOs.
            test(...
                [ 1, 3], [5, 8], {NSOID_1, NSOID_2}, ...
                [0:9], ...
                {...
                    logical([0, 1, 1, 1, 1, 1, 0, 0, 0, 0]'), ...
                    logical([0, 0, 0, 1, 1, 1, 1, 1, 1, 0]') ...
                }, ...
                {NSOID_1, NSOID_2}, [1, 2])

        end



        function test_read_file_BICAS(testCase)
            % NOTE: Only read BICAS's own default file (in BICAS's git repo).
            bicasRootPath = bicas.utils.get_BICAS_root_path();

            SETTINGS = bicas.create_default_SETTINGS();
            SETTINGS.make_read_only()
            rcsNsoRelativePath = SETTINGS.get_fv('PROCESSING.RCS_NSO.FILE.RELATIVE_PATH');

            nsoFilePath = fullfile(bicasRootPath, rcsNsoRelativePath);

            % TEST
            NsoTable = bicas.NSO_table.read_file_BICAS(nsoFilePath);
            testCase.verifyTrue(isa(NsoTable, 'bicas.NSO_table'))

            nEvents = irf.assert.sizes( ...
                NsoTable.evtStartTt2000Array, [-1, 1], ...
                NsoTable.evtStopTt2000Array,  [-1, 1], ...
                NsoTable.evtNsoIdCa,          [-1, 1]);
            testCase.verifyTrue(nEvents > 300)
        end



    end    % methods(Test)



end
