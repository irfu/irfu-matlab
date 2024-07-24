% Parse list of command-line options.
%
% Parses list of command-line arguments assuming it is composed of a list of
% options as defined below. The function tries to give accurate user-friendly
% errors (not assertions) for non-compliant arguments and absence of required
% arguments. The order of the options is arbitrary.
%
%
% IMPLEMENTATION NOTE: Reasons for using containers.Map (instead of arrays,
% array of structs, cell array, structure of structures).
% 1) Caller can easily build up list of options by amending list.
% 2) Can use key strings to identify options rather than the CLI option
% headers themselves (the latter could be short and cryptic, and can change
% since it is part of the command-line UI).
% 3) Can use more key strings (more permitted characters) than for structure
% field names.
% 4) Easy for caller to group subsets of options by keeping track of sets of
% keys. The caller can merge groups of options (before submitting as one
% parameter), and can split the returned result into the groups of options.
% Ex: Input files as opposed to output directory, log directory, config file.
% NOTE: The caller can easily(?) convert result into struct (one field per
% key-value pair).
%
%
% DEFINITIONS OF TERMS
% ====================
% Example argument list referred to below:
%       --verbose --file ~/bicas.conf --setting varX 123
% --
% Option header
%       A predefined (hard-coded, more or less) string meant to match a single
%       argument, e.g. "--verbose", "--file", "--setting". It does not have to
%       begin with "--" or "-" though that is the convention.
% Option value(s)
%       A sequence of arguments (could be zero arguments) following an option
%       header with which they are associated, e.g. none, "bicas.conf", or "X" &
%       "123". The number of expected option values should be predefined for the
%       option.
% Option
%       The combination of an option header and the immediately subsequent (and
%       associated) option values.
% Option ID
%       Unique, arbitrary string used to refer to the definition of an option
%       and the corresponding results from parsing CLI arguments.
%
%
% ARGUMENTS
% =========
% cliArgumentsCa
%       Column cell array of strings representing a sequence of CLI arguments.
% OptionsConfigMap
%       containers.Map
%       <keys>
%           Option ID.
%       <values>
%           One instance of class bicas.utils.cli.CliOptionConfig.
%
%
% RETURN VALUES
% =============
% OptionValuesMap
%       containers.Map with
%       <keys>
%         Option ID.
%       <values>
%         Column array of instances of class bicas.utils.cli.CliOptionValue.
%       NOTE: From this one can always read out whether an option was found
%       or not: even an option without option values contains a list of zero
%       values.
%       NOTE: Can read out the order of occurrence, e.g. for having a later
%       occurrence override a preceding one.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-02, reworked 2017-02-09, reworked 2019-07-22.
%
function OptionValuesMap = parse_CLI_options(cliArgumentsCa, OptionsConfigMap)
%
% PROPOSAL: Return some kind of help information to display proper user-friendly error message.
% PROPOSAL: Return struct array, one index for every option header+values (combined).
%   .index/.location : Number. Tells the order of (the groups of) arguments.
%   .optionId      :
%   .optionHeader  : String
%   .optionValues  : Cell array of strings
%   NOTE: Might want to sort/search by either .index or .optionId .
%   Therefore not Map. When there are several occurrences, then one can use
%   e.g. syntax
%           s=struct('x', {1,3,2,3}, 'y', {9,8,7,6})
%           xa = [s.x]; s(find(xa==3, 1, 'last')).y
%   CON: Return format is harder to search through when searching.
%       CON: Mostly/only when there are many occurrences.



% ASSERTIONS: Check argument types, sizes.
assert(iscell(cliArgumentsCa), 'cliArgumentsCa is not a cell array.')
assert(iscolumn(cliArgumentsCa))
assert(isa(OptionsConfigMap, 'containers.Map'))



[OptionsConfigMap, OptionValuesMap] = init_assert(OptionsConfigMap);
% Convert to struct array (NOT cell array of structs).
OptionsConfigArray = cellfun(@(x) (x), OptionsConfigMap.values);



%====================================
% Iterate over list of CLI arguments
%====================================
iCliArg = 1;
while iCliArg <= length(cliArgumentsCa)
  [OptionValuesMap, iCliArgLastValue] = try_interpret_option(...
    cliArgumentsCa, iCliArg, OptionsConfigArray, OptionValuesMap);

  iCliArg = iCliArgLastValue + 1;
end   % while



%=====================================================
% ASSERTION: Check that all required options were set
%=====================================================
for iOption = 1:length(OptionsConfigArray)
  optionId     = OptionsConfigArray(iOption).optionId;
  OptionConfig = OptionsConfigMap(optionId);
  optionValues = OptionValuesMap(optionId);

  if strcmp(OptionConfig.occurrenceRequirement, '0-1')

    if numel(optionValues) > 1
      error('BICAS:CLISyntax', ...
        'Found more than one occurrence of command-line option "%s".', ...
        OptionConfig.optionHeaderRegexp)
    end

  elseif strcmp(OptionConfig.occurrenceRequirement, '1')

    if numel(optionValues) ~= 1
      error('BICAS:CLISyntax', ...
        ['Could not find required command-line option matching', ...
        ' regular expression "%s".'], ...
        OptionConfig.optionHeaderRegexp)
    end

  elseif strcmp(OptionConfig.occurrenceRequirement, '0-inf')
    % Do nothing.

  else
    error('BICAS:Assertion', ...
      'Can not interpret OptionConfig.occurrenceRequirement="%s".', ...
      OptionConfig.occurrenceRequirement)

  end
end

end







% Try interpret a specific argument as an option header, followed by the
% expected number of option values.
%
% IMPLEMENTATION NOTE: Implemented as separate function to insulate the use of
% variables.
%
function [OptionValuesMap, iCliArgLastValue] = try_interpret_option(...
  cliArgumentsCa, iCliArg, OptionsConfigArray, OptionValuesMap)

cliArgument = cliArgumentsCa{iCliArg};

%=========================================
% Search for a matching CLI option string
%=========================================
% NOTE: It is more convenient to work with arrays than maps here.
iRegexpMatches  = find(irf.str.regexpf(...
  cliArgument, {OptionsConfigArray.optionHeaderRegexp}));

% IP = Interpretation Priority
% Array over regexp matches.
ipArray = [OptionsConfigArray(iRegexpMatches).interprPriority];
ip      = max(ipArray);    % Scalar
iMatch  = iRegexpMatches(ip == ipArray);



% ASSERTION
nMatchingOptions = numel(iMatch);
if nMatchingOptions == 0
  % CASE: Argument list does not conform to configuration.

  % NOTE: Phrase chosen for case that there may be multiple sequences of
  % arguments which are parsed separately.
  error('BICAS:CLISyntax', ...
    ['Can not interpret command-line argument "%s". It is not a', ...
    ' permitted option header in this sequence of arguments.'], ...
    cliArgument)
elseif nMatchingOptions >= 2
  % CASE: Configuration is bad (ill-defined).
  error('BICAS:Assertion', ...
    ['Can interpret CLI option in multiple ways, because the', ...
    ' interpretation of CLI arguments is badly configured.'])
end



%===========================================
% CASE: There is exacly one matching option
%===========================================
optionId = OptionsConfigArray(iMatch).optionId;

% Store values for this option. (May be zero option values).
OptionConfig = OptionsConfigArray(iMatch);
OptionValues = OptionValuesMap(optionId);

iCliArgLastValue = iCliArg + OptionConfig.nValues;
% ASSERTION: Argument list does not conform to configuration.
if iCliArgLastValue > length(cliArgumentsCa)
  error('BICAS:CLISyntax', ...
    ['Can not find the argument(s) that is/are expected to follow', ...
    ' command-line option header "%s".'], ...
    cliArgument)
end

% Extract option values associated with the option header.
OptionValues(end+1) = bicas.utils.cli.CliOptionValue(...
  iCliArg, ...
  cliArgumentsCa{iCliArg}, ...
  cliArgumentsCa(iCliArg+1:iCliArgLastValue, 1));

OptionValuesMap(optionId) = OptionValues;
end







% Various initializations and assertions.
%
% RETURN VALUES
% =============
% OptionsConfigMapModifCopy
%       Normalized and modified (added .optionId) deep copy of OptionsConfigMap.
% EmptyOptionValuesMap
%       containers.Map with
%       <Keys>
%         Same as OptionsConfigMap
%       <Values>
%         Empty struct array with pre-defined fields.
%         NOTE: Applies to both options with and without values!
%
function [OptionsConfigMapModifCopy, EmptyOptionValuesMap] = init_assert(...
  OptionsConfigMap)

% Copy filled with modified OptionsValuesMap.
EmptyOptionValuesMap      = containers.Map;
OptionsConfigMapModifCopy = containers.Map;

% List to iterate over map.
optionIdCa                = OptionsConfigMap.keys;
for iOption = 1:length(optionIdCa)
  optionId = optionIdCa{iOption};

  %===============================
  % Set OptionsConfigMapModifCopy
  %===============================
  assert(isa(OptionsConfigMap(optionId), 'bicas.utils.cli.CliOptionConfig'))
  ModifOptionConfig = struct(OptionsConfigMap(optionId));

  % ASSERTION: OptionConfig is the right struct.
  irf.assert.struct(ModifOptionConfig, ...
    {'optionHeaderRegexp', 'occurrenceRequirement', 'nValues'}, ...
    {'interprPriority'})

  % ASSERTION
  assert(isfinite(ModifOptionConfig.interprPriority))

  % Add .optionId
  ModifOptionConfig.optionId          = optionId;

  OptionsConfigMapModifCopy(optionId) = ModifOptionConfig;

  %==========================
  % Set EmptyOptionValuesMap
  %==========================
  EmptyOptionValuesMap(optionId) = bicas.utils.cli.CliOptionValue.empty(0, 1);
end

end
