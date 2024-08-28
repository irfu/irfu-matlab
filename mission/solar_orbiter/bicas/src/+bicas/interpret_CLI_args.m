%
% Function for doing the first interpretation of BICAS' CLI arguments and
% returns the data as a more easy-to-understand struct.
%
%
% ARGUMENTS
% =========
% cliArgumentsCa
%       Column cell array with all BICAS CLI arguments, both official and
%       unofficial.
%
%
% RETURN VALUE
% ============
% CliData
%       Struct with fields:
%       .bfmid               : String constant
%       .swmArg              : String constant
%       .icdLogFile          : Empty if argument not given.
%       .matlabLogFile       : Empty if argument not given.
%       .configFile          : Empty if argument not given.
%       .SipMap              : containers.Map with SIPs.
%                              key   = CLI argument without prefix
%                              value = file path (argument)
%       .ModifiedSettingsMap : containers.Map.
%                              key   = settings key (argument)
%                              value = settings value (argument)
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
% Reasons for having this function as separate function:
% -- Enable separate manual & automatic testing.
% -- Separate BICAS' "functionality" from "CLI syntax".
%    ==> Easier to change CLI syntax.
% -- Reduce size of BICAS main function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-22.
%
function CliData = interpret_CLI_args(cliArgumentsCa)
% PROPOSAL: Generic utility function for converting list of mutually
%       exclusive (assertion) booleans into one unique value.
%   Ex: Convert list of booleans for various argument flags (any
%       application) into one variable value.
%       Ex: Flags for BFMs.
%
% PROPOSAL: Use class for return value.

SWM_CLI_OPTION_REGEX = bicas.const.SWM_CLI_OPTION_REGEX;

%===============================================================================
% Configure
% (1) permitted RCS ICD CLI options COMMON for all BFMs
% (2) RCS ICD CLI options for SIPs
% (2) unofficial options
% NOTE: Exclude the SWM argument.
%===============================================================================
COPC_ARRAY = [...
  bicas.utils.cli.OptionConfig('VERSION_OPTION_ID',           '--version',          '0-1',   0,  0); ...
  bicas.utils.cli.OptionConfig('IDENTIFICATION_OPTION_ID',    '--identification',   '0-1',   0,  0); ...
  %bicas.utils.cli.OptionConfig('SWD_OPTION_ID',               '--swdescriptor',     '0-1',   0,  0); ...
  bicas.utils.cli.OptionConfig('HELP_OPTION_ID',              '--help',             '0-1',   0,  0); ...
  bicas.utils.cli.OptionConfig('SWM_OPTION_ID',               SWM_CLI_OPTION_REGEX, '0-1',   0,  0); ...
  ...
  % NOTE: ICD_LOG_FILE_OPTION_ID is an option to permit but ignore since it is handled by the bash launcher script, not the MATLAB code.
  bicas.utils.cli.OptionConfig('SIP_OPTION_ID',               '--(..*)',            '0-inf', 1, -1); ...
  bicas.utils.cli.OptionConfig('ICD_LOG_FILE_OPTION_ID',      '--log',              '0-1',   1,  0); ...
  bicas.utils.cli.OptionConfig('MATLAB_LOG_FILE_OPTION_ID',   '--log-matlab',       '0-1',   1,  0); ...
  bicas.utils.cli.OptionConfig('CONFIG_FILE_OPTION_ID',       '--config',           '0-1',   1,  0); ...
  ...
  % Unofficial arguments
  bicas.utils.cli.OptionConfig('MODIFIED_SETTINGS_OPTION_ID', '--set',              '0-inf', 2,  0); ...
  ];



CliData = [];



%============================================================================
% Extract the modified BICAS settings from the unofficial CLI arguments
% ---------------------------------------------------------------------
% IMPLEMENTATION NOTE: CliSettingsVsMap corresponds to one definition of ONE
% option (in the meaning of bicas.utils.cli.parse_CLI_options) and is filled
% with the corresponding option values in the order of the CLI arguments.
%   ==> A later occurrence of an option with the same first option
%       value, overwrites previous occurrences of the option with the same
%       first option value. This is the intended behaviour (not a side
%       effect).
%       Ex: --set SETTING_NAME 0 --setting SETTING_NAME 1
%============================================================================
CovcMap = bicas.utils.cli.parse_CLI_options(...
  cliArgumentsCa, COPC_ARRAY);
CliData.ModifiedSettingsMap = convert_modif_settings_COPVs_to_SettingsMap(...
  CovcMap('MODIFIED_SETTINGS_OPTION_ID'));



%=============================================================================
% Parse RCS ICD arguments
% -----------------------
% NOTE: Interprets RCS ICD as permitting (official) arguments next to non-SWM
% BFM mode arguments.
%=============================================================================

% Convert presence of BFM flag (mutually exclusive) into the correct constant.
% {i, 1} = false/true
% {i, 2} = BFMID
LogicalBfmidTable = {
  ~isempty(CovcMap('VERSION_OPTION_ID')),        'VERSION_BFM'; ...
  ~isempty(CovcMap('IDENTIFICATION_OPTION_ID')), 'IDENTIFICATION_BFM'; ...
  %~isempty(CovcMap('SWD_OPTION_ID')),            'SWD_BFM'; ...
  ~isempty(CovcMap('HELP_OPTION_ID')),           'HELP_BFM'; ...
  ~isempty(CovcMap('SWM_OPTION_ID')),            'SWM_BFM'};
assert(...
  sum([LogicalBfmidTable{:,1}]) == 1, ...
  'BICAS:interpret_CLI_syntax:CLISyntax', ...
  'Illegal combination of arguments.')
CliData.bfmid = LogicalBfmidTable{[LogicalBfmidTable{:, 1}], 2};

SipCovpArray = CovcMap('SIP_OPTION_ID');

switch CliData.bfmid

  %case {'VERSION_BFM', 'IDENTIFICATION_BFM', 'SWD_BFM', 'HELP_BFM'}
  case {'VERSION_BFM', 'IDENTIFICATION_BFM', 'HELP_BFM'}

    CliData.swmArg = [];
    CliData.SipMap = irf.ds.create_containers_Map('char', 'char', {}, {});
    assert(...
      isempty(SipCovpArray), ...
      'Specified illegal specific input parameters (SIP).')

  case 'SWM_BFM'

    CopvArray = CovcMap('SWM_OPTION_ID');

    % ASSERTION
    % NOTE: Somewhat of a hack, since can not read out from using
    % bicas.utils.cli.parse_CLI_options where the SWM option is located
    % among the arguments. The code knows it should be somewhere.
    if numel(CopvArray) ~= 1
      % Somewhat misleading error message. Hard to be accurate without
      % too much effort or by explaining the argument-parsing
      % algorithm to the user.
      error('BICAS:CLISyntax', 'Can not interpret argument(s).')
    elseif CopvArray.iOptionHeaderCliArgument ~= 1
      error('BICAS:CLISyntax', ...
        ['First argument can not be interpreted as', ...
        ' a S/W mode as expected.'])
    end

    CliData.swmArg = CopvArray.optionHeader;
    CliData.SipMap = convert_SIP_COPVs_to_Map(SipCovpArray);

  otherwise
    error('BICAS:Assertion', 'Illegal CliData.bfmid value.')
end



CovcArray = CovcMap('ICD_LOG_FILE_OPTION_ID');
if isempty(CovcArray)   CliData.icdLogFile = [];
else                    CliData.icdLogFile = CovcArray(end).optionValuesCa{1};
end

CovcArray = CovcMap('MATLAB_LOG_FILE_OPTION_ID');
if isempty(CovcArray)   CliData.matlabLogFile = [];
else                    CliData.matlabLogFile = CovcArray(end).optionValuesCa{1};
end

CovcArray = CovcMap('CONFIG_FILE_OPTION_ID');
if isempty(CovcArray)   CliData.configFile = [];
else                    CliData.configFile = CovcArray(end).optionValuesCa{1};
end



irf.assert.struct(CliData, ...
  {'bfmid', 'swmArg', 'icdLogFile', 'matlabLogFile', ...
  'configFile', 'SipMap', ...
  'ModifiedSettingsMap'}, {})

end



% NOTE: Checks (assertion) for doubles.
function SipMap = convert_SIP_COPVs_to_Map(CopvArray)
SipMap = irf.ds.create_containers_Map('char', 'char', {}, {});

for iSip = 1:numel(CopvArray)
  key = CopvArray(iSip).optionHeader(3:end);
  if SipMap.isKey(key)
    error('BICAS:CLISyntax', ...
      ['Specifying same specific input parameter (argument)', ...
      ' more than once.'])
  end
  SipMap(key) = CopvArray(iSip).optionValuesCa{1};
end
end



% NOTE: Deliberately does not check for doubles.
function SettingsMap = convert_modif_settings_COPVs_to_SettingsMap(CovcArray)
SettingsMap = irf.ds.create_containers_Map('char', 'char', {}, {});

for iSetting = 1:length(CovcArray)
  settingKey   = CovcArray(iSetting).optionValuesCa{1};
  settingValue = CovcArray(iSetting).optionValuesCa{2};
  SettingsMap(settingKey) = settingValue;
end
end
