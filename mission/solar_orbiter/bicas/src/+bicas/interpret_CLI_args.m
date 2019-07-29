%
% Function for doing the first interpretation of BICAS' CLI arguments and returns the data as a more easy-to-understand
% struct. This function covers both official and inofficial arguments.
%
%
% RETURN VALUE
% ============
% CliData : struct with fields:
%   .functionalityMode          : String constant
%   .swModeArg                  : String constant
%   .logFile                    : Empty if argument not given.
%   .configFile                 : Empty if argument not given.
%   .SpecInputParametersMap     : containers.Map. key=CLI argument without prefix --> value=file path
%   .ModifiedSettingsMap        : containers.Map. key=settings key                --> value=settings value
%
%
% IMPLEMENTATION NOTES
% ====================
% It is difficult to interpret BICAS arguments all in one go since some arguments do, or might, influence which other
% arguments (at least what the RCS ICD calls "specific inut parameters") are legal:
% -- functionality mode
% -- s/w mode
% -- config file        (might implement enabling/disabling s/w modes for debugging)
% -- settings arguments (might implement enabling/disabling s/w modes for debugging)
% To keep the function agnostic about s/w modes and input and output data sets, and settings, this function does not
% determine whether specific input parameters or settings keys/values) are legal. The caller has to do those checks.
%
%
% RATIONALE
% =========
% Reasons for having as separate function:
% -- Enable separate manual & automatic testing.
% -- Separate BICAS' "functionality" from "CLI syntax". ==> Easier to change CLI syntax.
% -- Reduce size of BICAS main function.
%
%
% TERMINOLOGY
% ===========
% FM = Functionality mode : Whether BICAS is launched with --version, --help, --identification, or (some) s/w mode.
% S/W mode                : Which dataset processing is to be performed.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-22.
%
function CliData = interpret_CLI_args(cliArgumentList)



C = bicas.error_safe_constants();



%==================================================================================
% Configure permitted RCS ICD CLI options COMMON for all BICAS functionality modes
% NOTE: Exclude the argument for functionality mode itself.
%==================================================================================
SW_MODE_REGEXP = '[^-][^-].*';

ICD_OPTIONS_CONFIG_MAP = containers.Map();
% NOTE: log_file and config_file are both options to permit but ignore since they are handled by bash launcher script.
ICD_OPTIONS_CONFIG_MAP('version_FM')                = struct('optionHeaderRegexp', '--version',        'occurrenceRequirement', '0-1',   'nValues', 0);
ICD_OPTIONS_CONFIG_MAP('identification_FM')         = struct('optionHeaderRegexp', '--identification', 'occurrenceRequirement', '0-1',   'nValues', 0);
ICD_OPTIONS_CONFIG_MAP('help_FM')                   = struct('optionHeaderRegexp', '--help',           'occurrenceRequirement', '0-1',   'nValues', 0);
ICD_OPTIONS_CONFIG_MAP('SW_mode')                   = struct('optionHeaderRegexp', SW_MODE_REGEXP,     'occurrenceRequirement', '0-1',   'nValues', 0);

ICD_OPTIONS_CONFIG_MAP('specific_input_parameters') = struct('optionHeaderRegexp', '--(.*)',           'occurrenceRequirement', '0-inf', 'nValues', 1, 'interprPriority', -1);

ICD_OPTIONS_CONFIG_MAP('log_file')                  = struct('optionHeaderRegexp', '--log',            'occurrenceRequirement', '0-1',   'nValues', 1);
ICD_OPTIONS_CONFIG_MAP('config_file')               = struct('optionHeaderRegexp', '--config',         'occurrenceRequirement', '0-1',   'nValues', 1);
% NOTE: "specific input parameter" is an RCS ICD term.

INOFF_OPTIONS_CONFIG_MAP = containers.Map();
INOFF_OPTIONS_CONFIG_MAP('modified_settings')       = struct('optionHeaderRegexp', '--set',    'occurrenceRequirement', '0-inf', 'nValues', 2);



%===============================================================================
% Separate CLI arguments into two different sequences/lists:
% (1) icdCliArgumentsList   = List of official arguments, as defined in RCS ICD.
% (2) inoffCliArgumentsList = List of inofficial arguments (may be empty).
%===============================================================================
iArgSeparator = find(strcmp(cliArgumentList, C.INOFFICIAL_ARGUMENTS_SEPARATOR));
if numel(iArgSeparator) == 0
    icdCliArgumentsList   = cliArgumentList;
    inoffCliArgumentsList = {};

elseif numel(iArgSeparator) == 1    % NOTE: Permit argument separator to be the very last argument.
    if (iArgSeparator <= 1)
        error('BICAS:CLISyntax', 'CLI argument separator at illegal position.')
    end
    icdCliArgumentsList   = cliArgumentList( 1 : (iArgSeparator-1) );
    inoffCliArgumentsList = cliArgumentList( (iArgSeparator+1) : end );
    
else
    error('BICAS:CLISyntax', 'Found more than one CLI argument separator.')
end



CliData = [];



%=======================================================================================================================
% Extract the modified settings from the inofficial CLI arguments.
%
% IMPLEMENTATION NOTE: CliSettingsVsMap corresponds to one definition of ONE option (in the meaning of
% parse_CLI_options) and is filled with the corresponding option values in the order of the CLI arguments.
%       ==> A later occurrence of an option with the same first option value, overwrites previous occurrences of the
%       option with the same first option value. This is the intended behaviour (not a side effect).
%       E.g. --setting LOGGING.IRF_LOG_LEVEL w --setting LOGGING.IRF_LOG_LEVEL n
%=======================================================================================================================
InoffOptionValuesMap = bicas.utils.parse_CLI_options(inoffCliArgumentsList, INOFF_OPTIONS_CONFIG_MAP);
CliData.ModifiedSettingsMap = convert_modif_settings_OptionValues_2_Map(InoffOptionValuesMap('modified_settings'));



%=====================================================================================================
% Parse RCS ICD arguments
% -----------------------
% NOTE: Interprets RCS ICD as permitting (official) arguments next to non-s/w mode functionality mode
% arguments.
%=====================================================================================================
CliData.SpecInputParametersMap = EJ_library.utils.create_containers_Map('char', 'char', {}, {});
IcdOptionValuesMap = bicas.utils.parse_CLI_options(icdCliArgumentsList, ICD_OPTIONS_CONFIG_MAP);


sipOptionValues = IcdOptionValuesMap('specific_input_parameters');
hasVersionFmOption = ~isempty(IcdOptionValuesMap('version_FM'));
hasIdentifFmOption = ~isempty(IcdOptionValuesMap('identification_FM'));
hasHelpFmOption    = ~isempty(IcdOptionValuesMap('help_FM'));
hasSwModeOption    = ~isempty(IcdOptionValuesMap('SW_mode'));

if     hasVersionFmOption && ~hasIdentifFmOption && ~hasHelpFmOption && ~hasSwModeOption
    CliData.functionalityMode = 'version';
    CliData.swModeArg         = [];
    assert(isempty(sipOptionValues), 'Specified illegal specific input parameters.')

elseif ~hasVersionFmOption && hasIdentifFmOption && ~hasHelpFmOption && ~hasSwModeOption
    CliData.functionalityMode = 'identification';
    CliData.swModeArg         = [];
    assert(isempty(sipOptionValues), 'Specified illegal specific input parameters.')

elseif ~hasVersionFmOption && ~hasIdentifFmOption && hasHelpFmOption && ~hasSwModeOption
    CliData.functionalityMode = 'help';
    CliData.swModeArg         = [];
    assert(isempty(sipOptionValues), 'Specified illegal specific input parameters.')

elseif ~hasVersionFmOption && ~hasIdentifFmOption && ~hasHelpFmOption && hasSwModeOption
    
    swModeArg = icdCliArgumentsList{1};    % NOTE: Always reading argument 1.
    
    % NOTE: Somewhat of a hack, since can not read out from using bicas.utils.parse_CLI_options where
    % the SW_mode option is located among the arguments. The code knows it should be somewhere.
    if ~EJ_library.utils.regexpf(swModeArg, SW_MODE_REGEXP)
        error('BICAS:CLISyntax', 'First argument can not be interpreted as a S/W mode as expected.')
    end

    CliData.functionalityMode = 'S/W mode';
    CliData.swModeArg         = swModeArg;   
    
    temp = convert_SIP_OptionValues_2_Map(sipOptionValues);
    CliData.SpecInputParametersMap = temp;
else
    error('BICAS:CLISyntax', 'Illegal combination of arguments.')
end



temp = IcdOptionValuesMap('log_file');
if isempty(temp)    CliData.logFile = [];
else                CliData.logFile = temp{end}{2};
end

temp = IcdOptionValuesMap('config_file');
if isempty(temp)    CliData.configFile = [];
else                CliData.configFile = temp{end}{2};
end

EJ_library.utils.assert.struct(CliData, {'functionalityMode', 'swModeArg', 'logFile', 'configFile', 'SpecInputParametersMap', 'ModifiedSettingsMap'})

end



% NOTE: Checks (assertion) for doubles.
function Map = convert_SIP_OptionValues_2_Map(optionValues)
Map = EJ_library.utils.create_containers_Map('char', 'char', {}, {});
for iSip = 1:numel(optionValues)
    temp = optionValues{iSip}{1};
    key = temp(3:end);
    if Map.isKey(key)
        error('BICAS:interpret_CLI_args:CLISyntax', 'Specifying same specific input parameter (argument) more than once.')
    end
    Map(key) = optionValues{iSip}{2};
end
end



% NOTE: Deliberately does not check for doubles.
function Map = convert_modif_settings_OptionValues_2_Map(optionValues)
Map = EJ_library.utils.create_containers_Map('char', 'char', {}, {});
for iSetting = 1:length(optionValues)
    Map(optionValues{iSetting}{2}) = optionValues{iSetting}{3};
end
end
