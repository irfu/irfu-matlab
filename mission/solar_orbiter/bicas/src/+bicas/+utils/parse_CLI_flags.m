% Parse list of command-line arguments (flags).
%
% Parses list of command-line arguments assuming it is a list of flags followed by a specified number of values.
% The function tries to give accurate user-friendly errors (not assertions) for non-compliant arguments and absence of
% required arguments.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02, reworked 2017-02-09
%
%
%
% DEFINITIONS OF TERMS
% ====================
% Flag   = (1) A predefined (hardcoded, more or less) string meant to match a single argument, e.g. "--verbose" or "--file".
%          (2) Flag (definition 1) plus its associated values.
% Values = A sequence of N arbitrary arguments that follow a flag (definition 1), e.g. "--file <value>".
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% cliArgumentsArray            : 1D cell array of strings.
% FlagsConfigMap               : containers.Map, number/string-->struct
%    <keys>                    : Arbitrary unique strings to identify the flags in the return value.
%    <values>                  : Struct. Information about each specified flag (syntax).
%       .cliFlagString             : The command-line flag string (e.g. "--version"), including any prefix (e.g. dash).
%       .occurranceRequirement : String specifying the number of times the flag may occurr. Alternatives: "0-1", "1", "0-inf".
%       .nValues               : The number of values that are expected after the flag.
%
% FlagValuesMap : containers.Map, number/string-->string/number
%    <keys>
%    <values> : Cell array of cell arrays of flag values, {iFlagOccurrance}{iValue}.
%               NOTE: From this one can always read out whether a flag was set: even a flag without values contains a list of zero values.
%               NOTE: A flag occurrence with zero values has an 1x0 cell array (not 0x0).
%
function FlagValuesMap = parse_CLI_flags(cliArgumentsArray, FlagsConfigMap)
%
% IMPLEMENTATION NOTE: Reasons for using containers.Map (instead of arrays, array of structs, cell
% array, structure of structures).
% 1) Caller can easily build up list of flags by amending list.
% 2) Can use key strings to identify flags rather than the CLI flag strings themselves (the latter could be short and
% cryptic).
% 3) Can use more key strings (more permitted characters) than for structure field names.
% 4) Easy for caller to group subsets of flags by keeping track of sets of keys. The caller can
% merge groups of flags (before submitting as one parameter), and can split the returned result into
% the groups of flags. Ex: Input files as opposed to output directory, log directory, config file.
% NOTE: The caller can easily(?) convert result into struct (one field per key-value pair).
%
% QUESTION: The function permits flags without value. Is this functionality really needed?



% ASSERTIONS
if ~iscell(cliArgumentsArray)
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a cell array.')
elseif length(cliArgumentsArray) ~= numel(cliArgumentsArray)
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a 1D cell array.')
elseif ~isa(FlagsConfigMap, 'containers.Map')
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a containers.Map.');
end



%===========================================
% Initializations before algorithm
%
% 1) Assertions
% 2) Collect cliFlagStringList
% 3) Initialize empty FlagValuesMap
%===========================================
cliFlagStringList = {};
FlagValuesMap = containers.Map;
flagKeysList    = FlagsConfigMap.keys;
flagConfigsList = FlagsConfigMap.values;
for iFlag = 1:length(flagKeysList)
    FlagConfig = flagConfigsList{iFlag};
    flagKey    = flagKeysList{iFlag};
    
    % ASSERTION
    if ~isempty( setxor(fieldnames(FlagConfig), {'cliFlagString', 'occurranceRequirement', 'nValues'}))
        error('parse_CLI_flags:Assertion:IllegalArgument', 'FlagsConfigMap does not have the expected fields.')
    end
    
    cliFlagStringList{end+1} = FlagConfig.cliFlagString;
    
    % Create return structure. Default: No flags found.
    % NOTE: Applies to both flags with and without values!
    FlagValuesMap(flagKey) = {};
end



% ASSERTION: Check that all cliFlagString are unique.
bicas.utils.assert_strings_unique(cliFlagStringList)



%====================================
% Iterate over list of CLI arguments
%====================================
iCliArgument = 1;
while iCliArgument <= length(cliArgumentsArray)
    cliArgument = cliArgumentsArray{iCliArgument};
    
    % Search for a matching CLI flag string. Set corresponding "flagKey".
    for iFlag = 1:length(flagKeysList)
        flagKey    = flagKeysList{iFlag};
        FlagConfig = flagConfigsList{iFlag};
        if strcmp(cliArgument, FlagConfig.cliFlagString)
            % CASE; Found CLI flag string.
            break
        else
            flagKey    = [];
            FlagConfig = [];   % Assignment is necessary in case loop ends without finding a matching a CLI flag key.
        end
    end
    if isempty(FlagConfig)
        error('parse_CLI_flags:CLISyntax', 'Can not interpret command-line argument "%s". There is no such flag.', cliArgument)
    end


    % Store values for this flag. (May be zero flag values).
    flagValues = FlagValuesMap(flagKey);
    iFirstValue = iCliArgument+1;
    iLastValue  = iCliArgument+FlagConfig.nValues;
    if iLastValue > length(cliArgumentsArray)
        error('parse_CLI_flags:CLISyntax', ...
            'Can not find the argument(s) that is/are expected to follow command-line flag "%s".', cliArgument)
    end
    flagValues{end+1} = cliArgumentsArray(iFirstValue:iLastValue);
    FlagValuesMap(flagKey) = flagValues;
    
    iCliArgument = iLastValue + 1;
end   % while



%========================================
% Check that all required flags were set
%========================================
for iFlag = 1:length(flagKeysList)
    flagKey    = flagKeysList{iFlag};
    FlagConfig = flagConfigsList{iFlag};
    flagValues = FlagValuesMap(flagKey);
    
    if strcmp(FlagConfig.occurranceRequirement, '0-1')
        if numel(flagValues) > 1
            error('parse_CLI_flags:CLISyntax', 'Found more than one occurance of command-line flag "%s".', FlagConfig.cliFlagString)
        end
    elseif strcmp(FlagConfig.occurranceRequirement, '1')
        if numel(flagValues) ~= 1
            error('parse_CLI_flags:CLISyntax', 'Could not find required command-line flag "%s".', FlagConfig.cliFlagString)
        end
    elseif strcmp(FlagConfig.occurranceRequirement, '0-N')
        ;   % Do nothing.
    else
        error('parse_CLI_flags:Assertion', 'Can not interpret occurranceRequirement="%s".', FlagConfig.occurranceRequirement)
    end
end

end
