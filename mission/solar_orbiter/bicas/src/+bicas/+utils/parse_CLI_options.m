% Parse list of command-line options.
%
% Parses list of command-line arguments assuming it is composed of a list of options as defined below. The function
% tries to give accurate user-friendly errors (not assertions) for non-compliant arguments and absence of required
% arguments. The order of the options is arbitrary.
%
%
% DEFINITIONS OF TERMS
% ====================
% Example argument list referred to below: --verbose --file ~/bicas.conf --setting varX 123
% --
% Option header   = A predefined (hardcoded, more or less) string meant to match a single argument, e.g. "--verbose",
%                   "--file", "--setting". It does not have to begin with "--" or "-" but that is the convention.
% Option value(s) = A sequence of arguments (could be zero arguments) following an option header with which they are
%                   associated, e.g. none, "bicas.conf", or "X" & "123". The number of expected option values should be
%                   predefined for the option.
% Option          = The combination of an option header and the subsequent option values.
% Option ID       = Unique, arbitrary string used to refer to the definition of an option and the corresponding results
%                   from parsing CLI arguments.
%
%
% ARGUMENTS
% =========
% cliArgumentsList             : 1D cell array of strings representing a sequence of CLI arguments.
% OptionsConfigMap             : containers.Map
%    <keys>                    : Option ID.
%    <values>                  : Struct. Information about each specified option (syntax).
%       .optionHeader          : The option header, including any prefix (e.g. dash).
%       .occurrenceRequirement : String specifying the number of times the option may occur.
%                                Permitted alternatives (strings):
%                                   '0-1'   = Option must occur once or never.
%                                   '1'     = Option must occur exactly once.
%                                   '0-inf' = Option may occur any number of times (zero or more).
%       .nValues               : The number of option values that must follow the option header.
%
%
% RETURN VALUES
% =============
% OptionValuesMap : containers.Map
%    <keys>   : Option ID.
%    <values> : Cell array of cell arrays of option values, {iOptionOccurrence}{iValue}.
%               NOTE: From this one can always read out whether an option was found or not: even an option without
%               option values contains a list of zero values.
%               NOTE: An option occurrence with zero values has a 1x0 cell array (not 0x0).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02, reworked 2017-02-09.
%
function OptionValuesMap = parse_CLI_options(cliArgumentsList, OptionsConfigMap)
%
% IMPLEMENTATION NOTE: Reasons for using containers.Map (instead of arrays, array of structs, cell
% array, structure of structures).
% 1) Caller can easily build up list of options by amending list.
% 2) Can use key strings to identify options rather than the CLI option headers themselves (the latter could be short
% and cryptic, and change).
% 3) Can use more key strings (more permitted characters) than for structure field names.
% 4) Easy for caller to group subsets of options by keeping track of sets of keys. The caller can
% merge groups of options (before submitting as one parameter), and can split the returned result into
% the groups of options. Ex: Input files as opposed to output directory, log directory, config file.
% NOTE: The caller can easily(?) convert result into struct (one field per key-value pair).



% ASSERTIONS: Check argument types, sizes.
if ~iscell(cliArgumentsList)
    error('parse_CLI_options:Assertion:IllegalArgument', 'Parameter is not a cell array.')
elseif length(cliArgumentsList) ~= numel(cliArgumentsList)
    error('parse_CLI_options:Assertion:IllegalArgument', 'Parameter is not a 1D cell array.')
elseif ~isa(OptionsConfigMap, 'containers.Map')
    error('parse_CLI_options:Assertion:IllegalArgument', 'Parameter is not a containers.Map.');
end



%===========================================
% Initializations before algorithm.
%
% 1) Assertions
% 2) Collect optionHeaderList
% 3) Initialize empty OptionValuesMap
%===========================================
optionHeaderList = {};
OptionValuesMap   = containers.Map;
optionIdsList     = OptionsConfigMap.keys;
optionConfigsList = OptionsConfigMap.values;
for iOption = 1:length(optionIdsList)
    OptionConfig = optionConfigsList{iOption};
    optionId     = optionIdsList{iOption};
    
    % ASSERTION: OptionConfig is the right struct.
    if ~isempty( setxor(fieldnames(OptionConfig), {'optionHeader', 'occurrenceRequirement', 'nValues'}))
        error('parse_CLI_options:Assertion:IllegalArgument', 'OptionsConfigMap does not have the expected fields.')
    end
    
    optionHeaderList{end+1} = OptionConfig.optionHeader;
    
    % Create return structure.
    % Default value: No options found.
    % NOTE: Applies to both options with and without values!
    OptionValuesMap(optionId) = {};
end



% ASSERTION: Check that all optionHeader are unique.
bicas.utils.assert_strings_unique(optionHeaderList)



%====================================
% Iterate over list of CLI arguments
%====================================
iCliArgument = 1;
while iCliArgument <= length(cliArgumentsList)
    cliArgument = cliArgumentsList{iCliArgument};
    
    % Search for a matching CLI option string. Set corresponding "optionId".
    optionId     = [];
    OptionConfig = [];
    for iOption = 1:length(optionIdsList)
        optionId = optionIdsList{iOption};
        if strcmp(cliArgument, optionConfigsList{iOption}.optionHeader)
            OptionConfig = optionConfigsList{iOption};
            % CASE: Found CLI option header.
            break
        end
    end
    if isempty(OptionConfig)
        % NOTE: Phrase chosen for case that there may be multiple sequences of arguments which are parsed separately.
        error('parse_CLI_options:CLISyntax', ...
            'Can not interpret command-line argument "%s". It is not a permitted option header in this sequence of arguments.', ...
            cliArgument)
    end


    % Store values for this option. (May be zero option values).
    optionValues = OptionValuesMap(optionId);
    iFirstValue = iCliArgument + 1;
    iLastValue  = iCliArgument + OptionConfig.nValues;
    if iLastValue > length(cliArgumentsList)
        error('parse_CLI_options:CLISyntax', ...
            'Can not find the argument(s) that is/are expected to follow command-line option header "%s".', cliArgument)
    end
    optionValues{end+1} = cliArgumentsList(iFirstValue:iLastValue);
    OptionValuesMap(optionId) = optionValues;
    
    iCliArgument = iLastValue + 1;
end   % while



%=====================================================
% ASSERTION: Check that all required options were set
%=====================================================
for iOption = 1:length(optionIdsList)
    optionId     = optionIdsList{iOption};
    OptionConfig = optionConfigsList{iOption};
    optionValues = OptionValuesMap(optionId);
    
    if strcmp(OptionConfig.occurrenceRequirement, '0-1')
        if numel(optionValues) > 1
            error('parse_CLI_options:CLISyntax', 'Found more than one occurance of command-line option "%s".', OptionConfig.optionHeader)
        end
    elseif strcmp(OptionConfig.occurrenceRequirement, '1')
        if numel(optionValues) ~= 1
            error('parse_CLI_options:CLISyntax', 'Could not find required command-line option "%s".', OptionConfig.optionHeader)
        end
    elseif strcmp(OptionConfig.occurrenceRequirement, '0-inf')
        ;   % Do nothing.
    else
        error('parse_CLI_options:Assertion', 'Can not interpret occurrenceRequirement="%s".', OptionConfig.occurrenceRequirement)
    end
end

end
