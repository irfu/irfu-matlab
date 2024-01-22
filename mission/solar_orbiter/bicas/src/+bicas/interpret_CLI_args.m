%
% Function for doing the first interpretation of BICAS' CLI arguments and
% returns the data as a more easy-to-understand struct. This function covers
% both official and unofficial arguments.
%
%
% RETURN VALUE
% ============
% CliData : struct with fields:
%   .functionalityMode          : String constant
%   .swmArg                  : String constant
%   .icdLogFile                 : Empty if argument not given.
%   .matlabLogFile              : Empty if argument not given.
%   .configFile                 : Empty if argument not given.
%   .SipMap                     : containers.Map with SIPs.
%                                   key   = CLI argument without prefix
%                                   value = file path (argument)
%   .ModifiedSettingsMap        : containers.Map.
%                                   key   = settings key (argument)
%                                   value = settings value (argument)
%
%
% IMPLEMENTATION NOTES
% ====================
% It is difficult to interpret BICAS arguments all in one go since some
% arguments do, or might, influence which other arguments (at least the SIPs)
% are legal:
% -- BICAS functionality mode (BFM)
% -- s/w mode (SWM)
% -- config file        (can enable/disable s/w modes)
% -- settings arguments (can enable/disable s/w modes)
% To keep the function agnostic about SWMs and input and output datasets, and
% settings, this function does not determine whether SIPs or settings
% keys/values) are legal. The caller has to do those checks.
%
%
% RATIONALE
% =========
% Reasons for having as separate function:
% -- Enable separate manual & automatic testing.
% -- Separate BICAS' "functionality" from "CLI syntax".
%    ==> Easier to change CLI syntax.
% -- Reduce size of BICAS main function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-22.
%
function CliData = interpret_CLI_args(cliArgumentList)
% PROPOSAL: Generic utility function for converting list of mutually
%       exclusive (assertion) booleans into one unique value.
%   Ex: Convert list of booleans for various argument flags (any
%       application) into one variable value.
%       Ex: Flags for BICAS functionality modes.
%
% PROPOSAL: Include assertion for unique input and output dataset paths.
%   NOTE: Assertion is presently in execute_SWM.
%
% PROPOSAL: Use classes instead of structs.

SWM_CLI_OPTION_REGEX = bicas.const.SWM_CLI_OPTION_REGEX;

%==================================================================================
% Configure
% (1) permitted RCS ICD CLI options COMMON for all BICAS functionality modes
% (2) RCS ICD CLI options for SIPs
% (2) unofficial options
% NOTE: Exclude the argument for functionality mode itself.
%==================================================================================
OPTIONS_CONFIG_MAP = containers.Map();
OPTIONS_CONFIG_MAP('VERSION_OPTION_ID')           = struct('optionHeaderRegexp', '--version',          'occurrenceRequirement', '0-1',   'nValues', 0);
OPTIONS_CONFIG_MAP('IDENTIFICATION_OPTION_ID')    = struct('optionHeaderRegexp', '--identification',   'occurrenceRequirement', '0-1',   'nValues', 0);
OPTIONS_CONFIG_MAP('SWD_OPTION_ID')               = struct('optionHeaderRegexp', '--swdescriptor',     'occurrenceRequirement', '0-1',   'nValues', 0);
OPTIONS_CONFIG_MAP('HELP_OPTION_ID')              = struct('optionHeaderRegexp', '--help',             'occurrenceRequirement', '0-1',   'nValues', 0);
OPTIONS_CONFIG_MAP('SWM_OPTION_ID')               = struct('optionHeaderRegexp', SWM_CLI_OPTION_REGEX, 'occurrenceRequirement', '0-1',   'nValues', 0);

% NOTE: ICD_LOG_FILE_OPTION_ID is an option to permit but ignore since it is handled by the bash launcher script, not the MATLAB code.
OPTIONS_CONFIG_MAP('SIP_OPTION_ID')               = struct('optionHeaderRegexp', '--(..*)',      'occurrenceRequirement', '0-inf', 'nValues', 1, 'interprPriority', -1);
OPTIONS_CONFIG_MAP('ICD_LOG_FILE_OPTION_ID')      = struct('optionHeaderRegexp', '--log',        'occurrenceRequirement', '0-1',   'nValues', 1);
OPTIONS_CONFIG_MAP('MATLAB_LOG_FILE_OPTION_ID')   = struct('optionHeaderRegexp', '--log-matlab', 'occurrenceRequirement', '0-1',   'nValues', 1);
OPTIONS_CONFIG_MAP('CONFIG_FILE_OPTION_ID')       = struct('optionHeaderRegexp', '--config',     'occurrenceRequirement', '0-1',   'nValues', 1);

% Unofficial arguments
OPTIONS_CONFIG_MAP('MODIFIED_SETTINGS_OPTION_ID') = struct('optionHeaderRegexp', '--set',        'occurrenceRequirement', '0-inf', 'nValues', 2);



CliData = [];



%============================================================================
% Extract the modified settings from the unofficial CLI arguments
% ---------------------------------------------------------------
% IMPLEMENTATION NOTE: CliSettingsVsMap corresponds to one definition of ONE
% option (in the meaning of parse_CLI_options) and is filled with the
% corresponding option values in the order of the CLI arguments.
%   ==> A later occurrence of an option with the same first option
%       value, overwrites previous occurrences of the option with the same
%       first option value. This is the intended behaviour (not a side
%       effect).
%       Ex: --set SETTING_NAME 0 --setting SETTING_NAME 1
%============================================================================
OptionValuesMap = bicas.utils.parse_CLI_options(...
  cliArgumentList, OPTIONS_CONFIG_MAP);
CliData.ModifiedSettingsMap = convert_modif_settings_OptionValues_2_Map(...
  OptionValuesMap('MODIFIED_SETTINGS_OPTION_ID'));



%=====================================================================
% Parse RCS ICD arguments
% -----------------------
% NOTE: Interprets RCS ICD as permitting (official) arguments next to
% non-s/w mode functionality mode arguments.
%=====================================================================
CliData.SipMap = irf.ds.create_containers_Map(...
  'char', 'char', {}, {});



sipOptionValues = OptionValuesMap('SIP_OPTION_ID');



% Convert presence of functionality mode flag (mutually exclusive) into the
% correct constant.
% {i, 1} = false/true
% {i, 2} = functionality mode string constant
tempTable = {
  ~isempty(OptionValuesMap('VERSION_OPTION_ID')),        'VERSION'; ...
  ~isempty(OptionValuesMap('IDENTIFICATION_OPTION_ID')), 'IDENTIFICATION'; ...
  ~isempty(OptionValuesMap('SWD_OPTION_ID')),            'SW_DESCRIPTOR'; ...
  ~isempty(OptionValuesMap('HELP_OPTION_ID')),           'HELP'; ...
  ~isempty(OptionValuesMap('SWM_OPTION_ID')),            'SW_MODE'};
assert(...
  sum([tempTable{:,1}]) == 1, ...
  'BICAS:interpret_CLI_syntax:CLISyntax', ...
  'Illegal combination of arguments.')
CliData.functionalityMode = tempTable{[tempTable{:, 1}], 2};

switch CliData.functionalityMode

  case {'VERSION', 'IDENTIFICATION', 'SW_DESCRIPTOR', 'HELP'}

    CliData.swmArg = [];
    assert(...
      isempty(sipOptionValues), ...
      'Specified illegal specific input parameters.')

  case 'SW_MODE'

    OptionValues = OptionValuesMap('SWM_OPTION_ID');

    % ASSERTION
    % NOTE: Somewhat of a hack, since can not read out from using
    % bicas.utils.parse_CLI_options where the SWM option is located
    % among the arguments. The code knows it should be somewhere.
    if numel(OptionValues) ~= 1
      % Somewhat misleading error message. Hard to be accurate without
      % too much effort or by explaining the argument-parsing
      % algorithm to the user.
      error('BICAS:CLISyntax', 'Can not interpret argument(s).')
    elseif OptionValues.iOptionHeaderCliArgument ~= 1
      error('BICAS:CLISyntax', ...
        ['First argument can not be interpreted as', ...
        ' a S/W mode as expected.'])
    end

    CliData.swmArg = OptionValues.optionHeader;

    CliData.SipMap = convert_SIP_OptionValues_2_Map(...
      sipOptionValues);

  otherwise
    error('BICAS:Assertion', 'Illegal CliData.functionalityMode value.')
end



temp = OptionValuesMap('ICD_LOG_FILE_OPTION_ID');
if isempty(temp)   CliData.icdLogFile = [];
else               CliData.icdLogFile = temp(end).optionValues{1};
end

temp = OptionValuesMap('MATLAB_LOG_FILE_OPTION_ID');
if isempty(temp)   CliData.matlabLogFile = [];
else               CliData.matlabLogFile = temp(end).optionValues{1};
end

temp = OptionValuesMap('CONFIG_FILE_OPTION_ID');
if isempty(temp)   CliData.configFile = [];
else               CliData.configFile = temp(end).optionValues{1};
end

irf.assert.struct(CliData, ...
  {'functionalityMode', 'swmArg', 'icdLogFile', 'matlabLogFile', ...
  'configFile', 'SipMap', ...
  'ModifiedSettingsMap'}, {})

end



% NOTE: Checks (assertion) for doubles.
function Map = convert_SIP_OptionValues_2_Map(optionValues)
Map = irf.ds.create_containers_Map('char', 'char', {}, {});

for iSip = 1:numel(optionValues)
  %temp = optionValues{iSip}{1};
  key = optionValues(iSip).optionHeader(3:end);
  %key = temp(3:end);
  if Map.isKey(key)
    error('BICAS:CLISyntax', ...
      ['Specifying same specific input parameter (argument)', ...
      ' more than once.'])
  end
  Map(key) = optionValues(iSip).optionValues{1};
end
end



% NOTE: Deliberately does not check for doubles.
function Map = convert_modif_settings_OptionValues_2_Map(optionValues)
Map = irf.ds.create_containers_Map('char', 'char', {}, {});

for iSetting = 1:length(optionValues)
  settingKey   = optionValues(iSetting).optionValues{1};
  settingValue = optionValues(iSetting).optionValues{2};
  Map(settingKey) = settingValue;
end
end
