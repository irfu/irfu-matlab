% ParsedCliArgumentsMap = parse_CLI_flags(cliArgumentsArray, FlagsConfigMap)   Parse list of command-line arguments (flags).
%
% Parses list of command-line arguments assuming it is a list of flags and pairs of flag+value.
% Function tries to give accurate user-friendly errors (not assertions) for non-compliant arguments and absence of
% required arguments.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
%
%
% DEFINITIONS OF TERMS
% ====================
% Flag  = A predefined (hardcoded, more or less) string meant to match a single argument, e.g. "--verbose" or "--file".
% Value = An argument that specifies a more or less arbitrary value that comes after some flags, e.g. "--file <value>".
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% cliArgumentsArray      : 1D cell array of strings.
% FlagsConfigMap : containers.Map, number/string-->struct
%    <keys>   : Arbitrary unique strings to identify the flags in the return value.
%    <values> : Information about each specified flag (syntax).
%       .cliString    : The command-line flag string (e.g. "--version"), including any prefix (e.g. dash).
%       .isRequired   : Whether the flag (and any specified value) is required as opposed to optional.
%       .expectsValue : Whether the flag expects the following argument to be a value connected to the flag.
%
% ParsedCliArgumentsMap : containers.Map, number/string-->string/number
%    <keys>
%    <values>  : For flag without value: 0=not set, 1=set
%                For flag with value:
%                      ~ischar(..) ==> the flag was not set (no such argument). (In reality numeric zero.)
%                       ischar(..) ==> the flag was set (incl. empty string). The string is the value.
%                   NOTE: Can not use/implement isempty(..) as criterion since empty strings can be
%                   legitimate argument values. Note that isempty('')==0.
%
function ParsedCliArgumentsMap = parse_CLI_flags(cliArgumentsArray, FlagsConfigMap)
%
% IMPLEMENTATION NOTE: Reasons for using containers.Map (instead of arrays, array of structs, cell
% array, structure of structures).
% 1) Caller can easily build up list of flags by amending list.
% 2) Can use key strings to identify flags (rather than the CLI flag strings themselves)
% 3) Can use more key strings (more characters) than for structure field names.
% 4) Easy for caller to group subsets of flags by keeping track of sets of keys. The caller can
% merge groups of flags (before submitting as one parameter), and can split the returned result into
% the groups of flags.
% Ex: Input files as opposed to output directory, log directory, config file.
% NOTE: The caller can easily(?) convert result into struct (one field per key-value pair).
%
% PROPOSAL: Change name: parse_arguments? parse_CLI_arguments?
%
% QUESTION: The function permits flags without value. Is this functionality really needed?
% PROPOSAL: Return containers.Map with values which are structures.
%    s.value
%    s.isSet   % True implies "value" if the corresponding flag has a value.



% ASSERTIONS
if ~iscell(cliArgumentsArray)
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a cell array.')
elseif length(cliArgumentsArray) ~= numel(cliArgumentsArray)
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a 1D cell array.')
elseif ~isa(FlagsConfigMap, 'containers.Map')
    error('parse_CLI_flags:Assertion:IllegalArgument', 'Parameter is not a containers.Map.');
end



cliStringArray = {};
ParsedCliArgumentsMap = containers.Map;
for skey = FlagsConfigMap.keys
    FlagConfig = FlagsConfigMap(skey{1});
    
    % ASSERTION
    if ~all(isfield(FlagConfig, {'cliString', 'isRequired', 'expectsValue'}))
        error('parse_CLI_flags:Assertion:IllegalArgument', 'FlagsConfigMap lacks all the required fields.')
    end
    
    cliStringArray{end+1} = FlagConfig.cliString;
    
    % Create return structure. Default: No flags found.
    % NOTE: Applies to both flags with and without values!
    ParsedCliArgumentsMap(skey{1}) = 0;         
end



% ASSERTION: Check that all cliString are unique.
bicas.utils.assert_strings_unique(cliStringArray)




%====================================
% Iterate over list of CLI arguments
%====================================
iCliArgument = 1;
while iCliArgument <= length(cliArgumentsArray)
    cliArgument = cliArgumentsArray{iCliArgument};
    
    % Find matching flag (and key). Set "flagKey".
    for skey = FlagsConfigMap.keys
        flagKey = skey{1};
        
        FlagConfig = FlagsConfigMap(flagKey);
        if strcmp(cliArgument, FlagConfig.cliString)
            break
        else
            FlagConfig = [];   % Assignment is necessary in case loop ends without finding a value.
        end
    end
    if isempty(FlagConfig)
        error('parse_CLI_flags:CLISyntax', 'Can not interpret command-line argument "%s". There is no such flag.', cliArgument)
    end
    
    
    
    % Check that flag has not been set already and that it does not.
    % NOTE: Can not use .isKey since all keys have been set to initialize ParsedCliArgumentsMap (needs default values
    %       for flags which are not set).
    if ~(isnumeric(ParsedCliArgumentsMap(flagKey)) && (ParsedCliArgumentsMap(flagKey) == 0))
        error('parse_CLI_flags:CLISyntax', 'The command-line flag "%s" was specified (at least) twice.', cliArgument)
    end
    
    if FlagConfig.expectsValue
        if iCliArgument >= length(cliArgumentsArray)
            error('parse_CLI_flags:CLISyntax', ...
                'Can not find the argument that is expected to follow command-line flag "%s".', cliArgument)
        end
        iCliArgument = iCliArgument + 1;
        ParsedCliArgumentsMap(flagKey) = cliArgumentsArray{iCliArgument};
    else
        ParsedCliArgumentsMap(flagKey) = 1;   % Set to "1" (true) representing that the flag was set.
    end
    
    iCliArgument = iCliArgument + 1;
end   % while



% Check that all required flags were set.
for skey = FlagsConfigMap.keys
    flag = FlagsConfigMap(skey{1});
    flagResult = ParsedCliArgumentsMap(skey{1});
    if flag.isRequired && (isnumeric(flagResult) && (flagResult == 0))
        % NOT: Not assertion since it is meant as a error message to be displayed for users.
        error('parse_CLI_flags:CLISyntax', 'Could not find required command-line flag "%s".', flag.cliString)
    end
end



% Convert flags (array) to structure.
% -----------------------------------
% IMPLEMENTATION NOTE: Up until here the algorithm works with the function argument describing the
% flags, and amends the results to it. First here does it copy the data to a variable that will be
% returned. This function has been modified several times to modify the output data structure and it
% is best to keep this conversion separate so that it is easy to modify the output data structure of
% the function again. The algorithm is also simpler(?) when working with arrays.
% flag_results = [];
% for jf = 1:length(flags)    
%     %flag_results.(flag.return_field_name).is_set = flags(jf).is_set;
%     %flag_results.(flag.return_field_name).value  = flags(jf).value;
%     
%     flag = flags(jf);
%     flag_results.(flag.return_field_name) = rmfield(flag, 'return_field_name');
%     
% end

end
