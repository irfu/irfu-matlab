% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-23
function interpret_CLI_args___ATEST

%     Ecm = EJ_library.utils.create_containers_Map('char', 'char', {}, {});

    tl = {};



    tl{end+1} = new_test( {'--help'}, ...
        'help',    [], [], [], {{}, {}}, {{}, {}});
    tl{end+1} = new_test( {'--version'}, ...
        'version', [], [], [], {{}, {}}, {{}, {}});
    tl{end+1} = new_test( {'--identification'}, ...
        'identification', [], [], [], {{}, {}}, {{}, {}});
    tl{end+1} = new_test( {'--swdescriptor'}, ...
        'S/W descriptor', [], [], [], {{}, {}}, {{}, {}});

    tl{end+1} = new_test( {'--help', '--log', 'logfile'}, ...
        'help', [], [], 'logfile', {{},{}}, {{},{}});
    tl{end+1} = new_test( {'--log', 'logfile', '--help'}, ...
        'help', [], [], 'logfile', {{},{}}, {{},{}});
    
    tl{end+1} = new_test( {'--version', '--log', 'logfile', '--config', 'configfile'}, ...
        'version', [], 'configfile', 'logfile', {{},{}}, {{},{}});

    tl{end+1} = new_test_EXC( {'--version', '--help'}, ...
        'MException');
    tl{end+1} = new_test_EXC( {'swmode', '--help'}, ...
        'MException');
    tl{end+1} = new_test_EXC( {'--in', 'infile', '--out', 'outfile', '--help'}, ...
        'MException');
    tl{end+1} = new_test_EXC( {'--in', 'infile', '--out', 'outfile', 'swmode'}, ...
        'MException');
    
    tl{end+1} = new_test( {'swmode', '--in', 'infile', '--out', 'outfile'}, ...
        'S/W mode', 'swmode', [], [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});
    
    tl{end+1} = new_test_EXC( {'--in', 'infile', 'swmode', '--out', 'outfile'}, ...
        'MException');    % S/w mode in the wrong place.
    

    
    tl{end+1} = new_test( {'swmode', '--in', 'infile', '--config', 'configfile', '--out', 'outfile'}, ...
        'S/W mode', 'swmode', 'configfile', [], {{}, {}}, {{'in', 'out'}, {'infile', 'outfile'}});
    
    % S/w mode in the wrong place.
    tl{end+1} = new_test_EXC( {'--in', 'infile', 'swmode', '--config', 'configfile', '--out', 'outfile'}, ...
        'MException');
    
    % S/w mode in TWO places (one correct). (Tested due to how algorithm works.)
    tl{end+1} = new_test_EXC( {'swmode1', '--in', 'infile', 'swmode2', '--config', 'configfile', '--out', 'outfile'}, ...
        'MException');
    
    tl{end+1} = new_test( {'--version'}, ...
                  'version', [], [], [], {{}, {}}, {{}, {}});
    tl{end+1} = new_test( {'--version', '--set', 'A', 'a', '--set', 'B', 'b'}, ...
                  'version', [], [], [], {{'A', 'B'}, {'a', 'b'}}, {{}, {}});



    EJ_library.atest.run_tests(tl)
end



% NOTE: Does not work when expecting an exception.
function Test = new_test(cliArgList, functionalityMode, swModeArg, configFile, logFile, ModifiedSettingsMap, SpecInputParametersMap)
outputs = {struct(...
    'functionalityMode', functionalityMode, ...
    'swModeArg',         swModeArg, ...
    'configFile',        configFile, ...
    'logFile',           logFile, ...
    'ModifiedSettingsMap',    EJ_library.utils.create_containers_Map('char', 'char', ModifiedSettingsMap{1},    ModifiedSettingsMap{2}), ...
    'SpecInputParametersMap', EJ_library.utils.create_containers_Map('char', 'char', SpecInputParametersMap{1}, SpecInputParametersMap{2})) ...
    };

assert(numel(outputs) == 1)

Test = EJ_library.atest.CompareFuncResult(@bicas.interpret_CLI_args, {cliArgList}, outputs);
end



% Test that generates exception.
% NOTE: Does not need arguments for outputs.
function Test = new_test_EXC(cliArgList, exceptionType)

Test = EJ_library.atest.CompareFuncResult(@bicas.interpret_CLI_args, {cliArgList}, exceptionType);
end



function OptionOccurrence = oo(iOptionHeaderCliArgument, optionHeader, optionValues)
    assert(iscell(optionValues))
    if isempty(optionValues)
        optionValues = cell(1,0);
    end
    OptionOccurrence = struct('iOptionHeaderCliArgument', iOptionHeaderCliArgument, 'optionHeader', optionHeader, 'optionValues', {optionValues});
end