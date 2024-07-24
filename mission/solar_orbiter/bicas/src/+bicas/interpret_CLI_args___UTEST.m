%
% matlab.unittest automatic test code for bicas.interpret_CLI_args().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12
%
classdef interpret_CLI_args___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Cycle through different filenames, paths, s/w mode strings,
  %           SIP option names.
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
        'HELP_BFM',    [], [], [], [], {{}, {}}, {{}, {}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--version'}, ...
        'VERSION_BFM', [], [], [], [], {{}, {}}, {{}, {}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--identification'}, ...
        'IDENTIFICATION_BFM', [], [], [], [], {{}, {}}, {{}, {}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--swdescriptor'}, ...
        'SWD_BFM', [], [], [], [], {{}, {}}, {{}, {}});



      %=======================
      % Non-s/w modes + extra
      %=======================
      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--help', '--log', 'ICD_LOG_FILE_OPTION_ID'}, ...
        'HELP_BFM', [], [], 'ICD_LOG_FILE_OPTION_ID', [], {{},{}}, {{},{}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--log', 'ICD_LOG_FILE_OPTION_ID', '--help'}, ...
        'HELP_BFM', [], [], 'ICD_LOG_FILE_OPTION_ID', [], {{},{}}, {{},{}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--log-matlab', 'MATLAB_LOG_FILE_OPTION_ID', '--help'}, ...
        'HELP_BFM', [], [], [], 'MATLAB_LOG_FILE_OPTION_ID', {{},{}}, {{},{}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--version', '--log', 'logfile', '--config', 'configfile'}, ...
        'VERSION_BFM', [], 'configfile', 'logfile', [], {{},{}}, {{},{}});

      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'--version', '--set', 'A', 'a', '--set', 'B', 'b'}, ...
        'VERSION_BFM', [], [], [], [], {{'A', 'B'}, {'a', 'b'}}, {{}, {}});

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
        'SWM_BFM', 'SWM', [], [], [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});

      % + extra
      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'SWM', '--in', 'infile', '--config', 'configfile', '--out', 'outfile'}, ...
        'SWM_BFM', 'SWM', 'configfile', [], [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});

      % Multiple input/output file arguments.
      bicas.interpret_CLI_args___UTEST.test(testCase, ...
        {'SWM', '--in1', 'infile1', '--out1', 'outfile1', '--in2', 'infile2', '--out2', 'outfile2'}, ...
        'SWM_BFM', 'SWM', [], [], [], {{}, {}}, ...
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
        cliArgCa, bfm, swmArg, configFile, ...
        icdLogFile, matlabLogFile, ...
        ModifiedSettingsMap, SipMap)

      cliArgCa = cliArgCa(:);

      assert(numel(ModifiedSettingsMap) == 2)
      assert(numel(SipMap)              == 2)

      ExpCliData = struct(...
        'bfm',                 bfm, ...
        'swmArg',              swmArg, ...
        'configFile',          configFile, ...
        'icdLogFile',          icdLogFile, ...
        'matlabLogFile',       matlabLogFile, ...
        'ModifiedSettingsMap', irf.ds.create_containers_Map(...
        'char', 'char', ModifiedSettingsMap{1:2} ...
        ), ...
        'SipMap', irf.ds.create_containers_Map(...
        'char', 'char', SipMap{1:2}) ...
        );
      assert(numel(ExpCliData) == 1, 'Bug in test code (not test case).')

      ActCliData = bicas.interpret_CLI_args(cliArgCa);

      testCase.verifyEqual(...
        ActCliData, ...
        ExpCliData)
    end



    function test_EXC(testCase, cliArgList)
      testCase.verifyError(...
        @() bicas.interpret_CLI_args(cliArgList), ...
        ?MException)
    end



  end    % methods(Static, Access=private)



end
