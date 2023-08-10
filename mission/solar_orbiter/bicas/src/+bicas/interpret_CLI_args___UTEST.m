%
% matlab.unittest automatic test code for bicas.interpret_CLI_args().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12
%
classdef interpret_CLI_args___UTEST < matlab.unittest.TestCase
    % PROPOSAL: Cycle through different filenames, paths, s/w mode strings,
    %           specific input parameter option names.
    %   CON: Most combinations provide no extra value. Only need to change one
    %        from a default test at a time.


    
    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test0(testCase)
            
            %=============== 
            % Non-s/w modes
            %===============
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--help'}, ...
                'help',    [], [], [], [], {{}, {}}, {{}, {}});
                        
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--version'}, ...
                'version', [], [], [], [], {{}, {}}, {{}, {}});
            
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--identification'}, ...
                'identification', [], [], [], [], {{}, {}}, {{}, {}});
            
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--swdescriptor'}, ...
                'S/W descriptor', [], [], [], [], {{}, {}}, {{}, {}});

            
            
            %=======================
            % Non-s/w modes + extra
            %=======================
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--help', '--log', 'ICD_log_file'}, ...
                'help', [], [], 'ICD_log_file', [], {{},{}}, {{},{}});
            
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--log', 'ICD_log_file', '--help'}, ...
                'help', [], [], 'ICD_log_file', [], {{},{}}, {{},{}});
            
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--log-matlab', 'MATLAB_log_file', '--help'}, ...
                'help', [], [], [], 'MATLAB_log_file', {{},{}}, {{},{}});

            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--version', '--log', 'logfile', '--config', 'configfile'}, ...
                'version', [], 'configfile', 'logfile', [], {{},{}}, {{},{}});
            
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'--version', '--set', 'A', 'a', '--set', 'B', 'b'}, ...
                'version', [], [], [], [], {{'A', 'B'}, {'a', 'b'}}, {{}, {}});
                        
            % ERRORS
            
            % Illegal combinations of modes.
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'--version', '--help'})
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'SWM', '--help'});
            
            
            
            %=========
            % S/w mode
            %=========
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'SWM', '--in', 'infile', '--out', 'outfile'}, ...
                'S/W mode', 'SWM', [], [], [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});
            
            % + extra
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'SWM', '--in', 'infile', '--config', 'configfile', '--out', 'outfile'}, ...
                'S/W mode', 'SWM', 'configfile', [], [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});

            % Multiple input/output file arguments.
            bicas.interpret_CLI_args___UTEST.test(testCase, ...
                {'SWM', '--in1', 'infile1', '--out1', 'outfile1', '--in2', 'infile2', '--out2', 'outfile2'}, ...
                'S/W mode', 'SWM', [], [], [], {{}, {}}, ...
                {{'in1', 'out1', 'in2', 'out2'}, ...
                 {'infile1', 'outfile1', 'infile2', 'outfile2'}});

            % ERRORS
            
            % S/w mode without s/w mode argument.
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'--in', 'infile', '--out', 'outfile', '--help'});
            
            % S/w mode, but s/w mode in the wrong place.
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'--in', 'infile', '--out', 'outfile', 'SWM'});
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'--in', 'infile', 'SWM', '--out', 'outfile'});
             
            % S/w mode in the wrong place.
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'--in', 'infile', 'SWM', '--config', 'configfile', '--out', 'outfile'});

            % S/w mode in TWO places (one correct). (Tested due to how algorithm works.)
            bicas.interpret_CLI_args___UTEST.test_EXC(testCase, ...
                {'SWM1', '--in', 'infile', 'SWM2', '--config', 'configfile', '--out', 'outfile'});

        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % NOTE: Does not work when expecting an exception.
        %
        % NOTE: testCase.verifyEqual() (as well as isequal()) do not care about
        % the order of key-vapue pairs in containers.Map.
        %
        function test(testCase, ...
                cliArgList, functionalityMode, swmArg, configFile, ...
                icdLogFile, matlabLogFile, ...
                ModifiedSettingsMap, SpecInputParametersMap)
            
            assert(numel(ModifiedSettingsMap)    == 2)
            assert(numel(SpecInputParametersMap) == 2)
            
            expOutput = struct(...
                'functionalityMode',      functionalityMode, ...
                'swmArg',              swmArg, ...
                'configFile',             configFile, ...
                'icdLogFile',             icdLogFile, ...
                'matlabLogFile',          matlabLogFile, ...
                'ModifiedSettingsMap',    irf.ds.create_containers_Map(...
                    'char', 'char', ModifiedSettingsMap{1:2}), ...
                'SpecInputParametersMap', irf.ds.create_containers_Map(...
                    'char', 'char', SpecInputParametersMap{1:2}) ...
            );
            assert(numel(expOutput) == 1, 'Bug in test code (not test case).')
            
            actOutput = bicas.interpret_CLI_args(cliArgList);
            
            testCase.verifyEqual(...
                actOutput, ...
                expOutput)
        end
        
        
        
        function test_EXC(testCase, cliArgList)
            testCase.verifyError(...
                @() bicas.interpret_CLI_args(cliArgList), ...
                ?MException)
        end
        
        
        
    end    % methods(Static, Access=private)

    
    
end
