%
% matlab.unittest automatic test code for bicas.utils.parse_CLI_options().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-10
%
classdef parse_CLI_options___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            % NOTE: Arguments OptionsConfigMap, inputStr switch places to make
            %       test code look better.
            function test(OptionsConfigMap, inputStr, outputMapKeys, outputMapValues)
                cliArgumentsList = strsplit(inputStr);

                expOutput = containers.Map(outputMapKeys, outputMapValues);

                actOutput = bicas.utils.parse_CLI_options(...
                    cliArgumentsList, OptionsConfigMap);

                testCase.verifyEqual(actOutput, expOutput)
            end



            % NOTE: Arguments OptionsConfigMap, inputStr switch places to make
            %       test code look better.
            function test_exc(OptionsConfigMap, inputStr)
                cliArgumentsList = strsplit(inputStr);

                testCase.verifyError(...
                    @() bicas.utils.parse_CLI_options(...
                        cliArgumentsList, OptionsConfigMap), ...
                        ?MException)
            end

            %===================================================================
            import bicas.utils.parse_CLI_options___UTEST.ocme
            import bicas.utils.parse_CLI_options___UTEST.oo



            % OCM = Options Config Map
            OCM1 = containers.Map(...
                {'a', 'b', 'c'}, ...
                {...
                    ocme('-a', '0-1',   0), ...
                    ocme('-b', '1',     1), ...
                    ocme('-c', '0-inf', 2)...
                });
            OCM2 = containers.Map(...
                {'--', '=='}, ...
                {...
                    ocme('--.*', '0-1',   0), ...
                    ocme('==.*', '0-inf', 0) ...
                });
            OCM3 = containers.Map(...
                {'all', 'log', 'set'}, ...
                {...
                    ocme('--.*',    '0-inf', 1, -1), ...
                    ocme('--log',   '1',     1), ...
                    ocme('--set.*', '0-inf', 1) ...
                });

            EOO = oo(cell(0,1), cell(0,1), cell(0,1));



            % Missing option -b.
            test_exc(OCM1, '-a')



            test(OCM1, '-b 123',           {'a', 'b', 'c'}, {EOO,               oo(1, '-b', {'123'}),   EOO});
            test(OCM1, '-a -b 123',        {'a', 'b', 'c'}, {oo(1, '-a', {}),   oo(2, '-b', {'123'}),   EOO});
            test(OCM1, '-a -b 123 -c 8 9', {'a', 'b', 'c'}, {oo(1, '-a', {}),   oo(2, '-b', {'123'}),   oo(4, '-c', {'8', '9'})});



            % Test multiple occurrences of the same option.
            test(OCM1, '-c 6 7 -a -b 123 -c 8 9', ...
                {'a', 'b', 'c'}, ...
                { oo(4, '-a', {}), ...
                  oo(5, '-b', {'123'}), ...
                 [oo(1, '-c', {'6', '7'}), ...
                  oo(7, '-c', {'8', '9'})...
                 ]});
             % Test multiple occurrences of the same option.
            test(OCM1, '-c 6 7 -b 123 -c 8 9', ...
                {'a', 'b', 'c'}, ...
                {EOO, ...
                 oo(4, '-b', {'123'}), ...
                [oo(1, '-c', {'6', '7'}), ...
                 oo(6, '-c', {'8', '9'})...
                 ]});

            test(OCM2, '--ASD',           {'--', '=='}, ...
                {oo(1, '--ASD', {}), ...
                 EOO});
            test(OCM2, '==ASD ==a --abc', {'--', '=='}, ...
                { oo(3, '--abc', {}), ...
                 [oo(1, '==ASD', {}), ...
                  oo(2, '==a',   {})...
                  ]});



            test(OCM3, '--input1 i1 --output1 o1 --log logfile', ...
                {'all', 'log', 'set'}, ...
                {[oo(1, '--input1',  {'i1'}), ...
                  oo(3, '--output1', {'o1'})...
                  ], ...
                  oo(5, '--log',     {'logfile'}), ...
                  EOO});

            test(OCM3, '--input1 i1 --output1 o1 --log logfile --setDEBUG ON', ...
                {'all', 'log', 'set'}, ...
                {[oo(1, '--input1',  {'i1'}), ...
                  oo(3, '--output1', {'o1'})...
                  ], ...
                  oo(5, '--log',      {'logfile'}), ...
                  oo(7, '--setDEBUG', {'ON'})});

        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        function Ocme = ocme(optionHeaderRegexp, occurrenceRequirement, nValues, interprPriority)
            args = {...
                'optionHeaderRegexp',    optionHeaderRegexp, ...
                'occurrenceRequirement', occurrenceRequirement, ...
                'nValues',               nValues};
            if nargin == 4
                args(end+1:end+2) = {...
                    'interprPriority', interprPriority};
            end

            Ocme = struct(args{:});
        end



        function OptionOccurrence = oo(...
                iOptionHeaderCliArgument, optionHeader, optionValues)

            assert(iscell(optionValues))
            if isempty(optionValues)
                % CASE: Option was never used.
                optionValues = cell(0,1);
            end
            optionValues = optionValues(:);   % Force column vector.

            OptionOccurrence = struct(...
                'iOptionHeaderCliArgument', iOptionHeaderCliArgument, ...
                'optionHeader',             optionHeader, ...
                'optionValues',             {optionValues});
        end



    end    % methods(Static, Access=private)



end
