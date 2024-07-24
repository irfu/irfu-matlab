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
%
%
% ARGUMENTS
% =========
% cliArgumentsCa
%       Column cell array of strings representing a sequence of CLI arguments.
% CopcArray
%       Column array of instances of COPC.
%
%
% RETURN VALUES
% =============
% CopvMap
%       containers.Map with
%       <keys>
%         Option ID.
%       <values>
%         Column array of instances of COPV.
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
function CopvMap = parse_CLI_options(cliArgumentsCa, CopcArray)
%
% PROPOSAL: Return some kind of help information to display proper user-friendly error message.



% ASSERTIONS: Check argument types, sizes.
assert(iscell(cliArgumentsCa), 'cliArgumentsCa is not a cell array.')
assert(iscolumn(cliArgumentsCa))
assert(isa(CopcArray, 'bicas.utils.cli.OptionConfig') & iscolumn(CopcArray))
irf.assert.castring_set({CopcArray.optionId})



% Create CopvMap: containers.Map
CopvMap = containers.Map;
for iOption = 1:length(CopcArray)
  optionId = CopcArray(iOption).optionId;
  CopvMap(optionId) = bicas.utils.cli.OptionValue.empty(0, 1);
end



%====================================
% Iterate over list of CLI arguments
%====================================
iCliArg = 1;
while iCliArg <= length(cliArgumentsCa)
  % NOTE: MODIFIES THE ARGUMENT "CopvMap" (handle object)!
  iCliArgLastValue = try_interpret_option(...
    cliArgumentsCa, iCliArg, CopcArray, CopvMap);

  iCliArg = iCliArgLastValue + 1;
end   % while



%=====================================================
% ASSERTION: Check that all required options were set
%=====================================================
for iOption = 1:length(CopcArray)
  optionId  = CopcArray(iOption).optionId;
  Copc      = CopcArray(iOption);
  CopvArray = CopvMap(optionId);

  if strcmp(Copc.occurrenceRequirement, '0-1')

    if numel(CopvArray) > 1
      error('BICAS:CLISyntax', ...
        'Found more than one occurrence of command-line option "%s".', ...
        Copc.optionHeaderRegexp)
    end

  elseif strcmp(Copc.occurrenceRequirement, '1')

    if numel(CopvArray) ~= 1
      error('BICAS:CLISyntax', ...
        ['Could not find required command-line option matching', ...
        ' regular expression "%s".'], ...
        Copc.optionHeaderRegexp)
    end

  elseif strcmp(Copc.occurrenceRequirement, '0-inf')
    % Do nothing.

  else
    error('BICAS:Assertion', ...
      'Can not interpret OptionConfig.occurrenceRequirement="%s".', ...
      Copc.occurrenceRequirement)

  end
end

end







% Try interpret a specific argument as an option header, followed by the
% expected number of option values.
%
% IMPLEMENTATION NOTE: Implemented as separate function to insulate the use of
% variables.
%
% ARGUMENTS
% =========
% CopvMap
%     NOTE: MODIFIES THIS ARGUMENT (handle object).
%
function iCliArgLastValue = try_interpret_option(...
  cliArgumentsCa, iCliArg, CopcArray, CopvMap)

cliArgument = cliArgumentsCa{iCliArg};

%=========================================
% Search for a matching CLI option string
%=========================================
% NOTE: It is more convenient to work with arrays than containers.Map here.
iRegexpMatches  = find(irf.str.regexpf(...
  cliArgument, {CopcArray.optionHeaderRegexp}));

% IP = Interpretation Priority
% Array over regexp matches.
ipArray = [CopcArray(iRegexpMatches).interprPriority];
ip      = max(ipArray);    % Scalar
iMatch  = iRegexpMatches(ip == ipArray);



% ASSERTION: There is exactly one matching (after priority) COPC.
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
optionId = CopcArray(iMatch).optionId;

% Store values for this option. (May be zero option values).
Copc      = CopcArray(iMatch);
CopvArray = CopvMap(optionId);

iCliArgLastValue = iCliArg + Copc.nValues;
% ASSERTION: Argument list does not conform to configuration.
if iCliArgLastValue > length(cliArgumentsCa)
  error('BICAS:CLISyntax', ...
    ['Can not find the argument(s) that is/are expected to follow', ...
    ' command-line option header "%s".'], ...
    cliArgument)
end

% Extract option values associated with the option header.
CopvArray(end+1) = bicas.utils.cli.OptionValue(...
  iCliArg, ...
  cliArgumentsCa{iCliArg}, ...
  cliArgumentsCa(iCliArg+1:iCliArgLastValue, 1));

% NOTE: MODIFIES THIS ARGUMENT (handle object).
CopvMap(optionId) = CopvArray;
end
