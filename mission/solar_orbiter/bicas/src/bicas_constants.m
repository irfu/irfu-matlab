% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software.
%
function constants = bicas_constants()
%
% IMPORTANT NOTE: Some constants (1) correspond exactly to fields in the S/W (JSON) descriptor, and
% (2) are unlikely to be used for anything else. These are labeled with a prefix "SWD_". Other
% variables which do not have the prefix may also be used for the S/W descriptor too but they are
% probably more unambiguous in their meaning.
% 
% IMPLEMENTATION NOTE: There is other code which builds a structure corresponding to the S/W descriptor
% from the constants structure here.
% Reasons for not putting the S/W descriptor structure inside the constants structure:
% (1) Some of the S/W descriptor variables have vague or misleading names (name, dataset versions, dataset IDs)
% (2) Some of the S/W descriptor variables are grouped in a way which
%     does not fit the rest of the code (modes[].outputs.output_XX.release).
% (2) Some of the S/W descriptor values are really structure field NAMES, but should be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (3) The constants structure would become dependent on the format of the S/W descriptor structure.
%     The latter might change in the future.
% (4) Some S/W descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (5) It is easier to add automatic checks on the S/W descriptor in the code that derives it.
%
% PROPOSAL: Use ~dataset ID for sw_mode.inputs{end}.CLI_parameter_name?
% PROPOSAL: Add root path?! How?
%    CON: Not hardcoded.
%
% PROPOSAL: Merge with init global constants?
%    CON: Calling would slow down code. Many functions use ERROR_CODES.
%       PROPOSAL: Use one big global variable instead (initialized here).
%    CON: Would require to having only "error safe code".
%       PROPOSAL: Limit to code without branching (except lazy eval.) and without reading/calling
%                 external files so that one single clear+run will tell whether it works or not.
%          Ex: for loops are then safe.
%          CON: Can not use irf.log messages (for temp. data).
% PROPOSAL: Merge with get_sw_descriptor?
%    CON: Adds much code (non-hardcoded constants).
%       CON: Can separate into separate functions.
%          PROPOSAL: root function: lazy eval.+call get_constants+add SWD with get_sw_descriptor+call validate.
%                    get_constants : Only ~hardcoded constants
% PROPOSAL: Add S/W descriptor to this structure by calling get_sw_descriptor?
%    PRO: Can force validation taking place in get_sw_descriptor.
% 
% PROPOSAL: Could validate only when extra parameters (e.g. .S/W root directory) are given.
% PROPOSAL: Validation in separate function (other file). Called separately.
%    PRO: Validation may require extra parameters. Ex: S/W root directory
%    PRO: Error-prone validation code can then be inside main try-catch.
% PROPOSAL: More validation: Check that master cdfs exist, that executable exists.



% Support caching/lazy evaluation.
persistent C
if ~isempty(C)
    constants = C;
    return
end



INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';
TEST_MASTER_CDF_FILENAME = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E-XXXTESTXXX_V01.cdf';   % For testing.



C.author_name = 'Erik P G Johansson';
C.author_email = 'erik.johansson@irfu.se';
C.institute = 'IRF-U';
C.master_cdfs_dir_rel = 'data';    % Location of master CDF files. Relative to the software directory structure root.


irf.log('w', 'Using temporary S/W name in constants.')
C.SWD_identification.project     = 'ROC-SGSE';
C.SWD_identification.name        = 'BICAS (temporary name)';   % Temporary sw name
C.SWD_identification.identifier  = 'ROC-SGSE-BICAS';           % Temporary sw name
C.SWD_identification.description = 'BIAS Calibration Software (BICAS) which derives the BIAS L2S input signals (plus some) from the BIAS L2R output signals.';

% Refers to the S/W descriptor release data for the entire software (not specific outputs).
C.SWD_release.version      = '0.0.1';
C.SWD_release.date         = '2016-06-03';
C.SWD_release.author       = C.author_name;
C.SWD_release.contact      = C.author_email;
C.SWD_release.institute    = C.institute;
C.SWD_release.modification = INITIAL_RELEASE_MODIFICATION_STR;

C.SWD_environment.executable = 'roc/bicas';   % Temporary sw name

C.sw_modes = {};

%---------------------------------------------------------------------------------------------------
% NOTE: sw_mode.inputs{..}.CLI_parameter_name = CLI parameter without flag prefix ("--").

sw_mode = [];
%sw_mode.CLI_parameter = 'L2S_LFR-SURV-CWF-E';
sw_mode.CLI_parameter = 'LFR-CWF-E';
sw_mode.SWD_purpose = 'Generate continuous waveform electric field data (potential difference) from LFR';

sw_mode.inputs = {[]};
sw_mode.inputs{end}.CLI_parameter_name  = 'input_hk';
sw_mode.inputs{end}.dataset_ID          = 'ROC-SGSE_HK_RPW-BIA';
sw_mode.inputs{end}.dataset_version_str = '01';

sw_mode.inputs{end+1} = [];
sw_mode.inputs{end}.CLI_parameter_name  = 'input_sci';
sw_mode.inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF';
sw_mode.inputs{end}.dataset_version_str = '01';

sw_mode.outputs = {[]};
sw_mode.outputs{end}.JSON_output_file_identifier = 'output_SCI';
sw_mode.outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
sw_mode.outputs{end}.dataset_version_str         = '01';
%sw_mode.outputs{end}.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01.cdf';
sw_mode.outputs{end}.master_cdf_filename         = TEST_MASTER_CDF_FILENAME;
sw_mode.outputs{end}.SWD_name                    = 'LFR L2s CWF science electric data in survey mode';
sw_mode.outputs{end}.SWD_description             = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
sw_mode.outputs{end}.SWD_level                   = 'L2S';
sw_mode.outputs{end}.SWD_release_date            = '2016-06-03';
sw_mode.outputs{end}.SWD_release_modification    = INITIAL_RELEASE_MODIFICATION_STR;

C.sw_modes{end+1} = sw_mode;
%---------------------------------------------------------------------------------------------------
sw_mode = [];
%sw_mode.CLI_parameter = 'L2S_LFR-SURV-SWF-E';
sw_mode.CLI_parameter = 'LFR-SWF-E';
sw_mode.SWD_purpose = 'Generate snapshow waveform electric (potential difference) data from LFR';

sw_mode.inputs = {[]};
sw_mode.inputs{end}.CLI_parameter_name  = 'input_hk';
sw_mode.inputs{end}.dataset_ID          = 'ROC-SGSE_HK_RPW-BIA';
sw_mode.inputs{end}.dataset_version_str = '01';

sw_mode.inputs{end+1} = [];
sw_mode.inputs{end}.CLI_parameter_name  = 'input_sci';
sw_mode.inputs{end}.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF';
sw_mode.inputs{end}.dataset_version_str = '01';

sw_mode.outputs = {[]};
sw_mode.outputs{end}.JSON_output_file_identifier = 'output_SCI';
sw_mode.outputs{end}.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
sw_mode.outputs{end}.dataset_version_str         = '01';
%sw_mode.outputs{end}.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01.cdf';
sw_mode.outputs{end}.master_cdf_filename         = TEST_MASTER_CDF_FILENAME;
sw_mode.outputs{end}.SWD_name                    = 'LFR L2s SWF science electric data in survey mode';
sw_mode.outputs{end}.SWD_description             = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
sw_mode.outputs{end}.SWD_level                   = 'L2S';
sw_mode.outputs{end}.SWD_release_date            = '2016-06-03';
sw_mode.outputs{end}.SWD_release_modification    = INITIAL_RELEASE_MODIFICATION_STR;

C.sw_modes{end+1} = sw_mode;
%---------------------------------------------------------------------------------------------------

validate_constants(C)

constants = C;

end
%===================================================================================================
% Any code for double-checking the validity of hardcoded constants.
%
function validate_constants(C)

global ERROR_CODES

for sw_mode = C.sw_modes
    % The RCS ICD, iss2rev2, section 3.2.3 only permits these characters.
    perm_chars = 'abcdefghijklmnopqrstuvxyz0123456789_';
    
    for input = sw_mode{1}.inputs
        disallowed_chars = setdiff(input{1}.CLI_parameter_name, perm_chars);
        if ~isempty(disallowed_chars)
            errorp(ERROR_CODES.ASSERTION_ERROR, 'Constants value contains illegal characters. This indicates a pure bug.');
        end
    end
end

end