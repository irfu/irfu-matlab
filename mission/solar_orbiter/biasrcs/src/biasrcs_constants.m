% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Defines constants used by the software.
%
function constants = biasrcs_constants()
%
% IMPORTANT NOTE: Some constants (1) correspond exactly to fields in the SW (JSON) descriptor,
% and (2) are unlikely to be used for anything else. These are labeled with a prefix
% "SWD_".
% 
% IMPLEMENTATION NOTE: There is other code which builds a structure corresponding to the SW descriptor
% from the constants structure here.
% Reasons for not putting the SW descriptor structure inside the constants structure:
% (1) Some of the SW descriptor variables have vague or misleading names (name, dataset versions, dataset IDs)
% (2) Some of the SW descriptor variables are grouped in a way which
%     does not fit the rest of the code (modes[].outputs.output_XX.release).
% (2) Some of the SW descriptor values are really structure field NAMES, but should be better as
%     structure field VALUES (e.g. input CLI parameter, output JSON identifier string).
% (3) The constants structure would become dependent on the format of the SW descriptor structure.
%     The latter might change in the future.
% (4) Some SW descriptor values repeat or can be derived from other constants (e.g. author, contact,
%     institute, output_XX.release.file).
% (5) It is easier to add automatic checks on the SW descriptor in the code that derives it.

persistent C
if ~isempty(C)
    constants = C;
    return
end

INITIAL_RELEASE_MODIFICATION_STR = 'No modification (initial release)';



C.author_name = 'Erik P G Johansson';
C.author_email = 'erik.johansson@irfu.se';
C.institute = 'IRF-U';



C.SWD_identification.project     = 'ROC-SGSE';
C.SWD_identification.name        = 'BIASRCS (temporary name)';   % Temporary sw name
C.SWD_identification.identifier  = 'ROC-SGSE-BIASRCS';           % Temporary sw name
C.SWD_identification.description = 'Calibration software which derives the BIAS input signals (plus some) from the BIAS output signals.';

% Refers to the S/W descriptor release data for the entire software (not specific outputs).
C.SWD_release.version      = '0.0.1';
C.SWD_release.date         = '2016-06-01';
C.SWD_release.author       = C.author_name;
C.SWD_release.contact      = C.author_email;
C.SWD_release.institute    = C.institute;
C.SWD_release.modification = INITIAL_RELEASE_MODIFICATION_STR;

C.SWD_environment.executable = 'roc/biasrcs';   % Temporary sw name

C.sw_modes = {};

%---------------------------------------------------------------------------------------------------
sw_mode = [];
%sw_mode.CLI_parameter = 'L2S_LFR-SURV-CWF-E';
sw_mode.CLI_parameter = 'LFR-CWF-E';
sw_mode.SWD_purpose = 'Generate continuous waveform electric field data (potential difference) from LFR';

sw_mode.inputs = {};

input = [];
input.CLI_parameter_name  = 'input_HK';
input.dataset_ID          = 'ROC-SGSE_HK_RPW-BIA';
input.dataset_version_str = '01';
sw_mode.inputs{end+1} = input;

input = [];
input.CLI_parameter_name  = 'input_SCI';
input.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF';
input.dataset_version_str = '01';
sw_mode.inputs{end+1} = input;

sw_mode.outputs = {};

output = [];
output.JSON_output_file_identifier = 'output_SCI';
output.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E';
output.dataset_version_str         = '01';
output.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01.cdf';
output.SWD_name                 = 'LFR L2s CWF science electric data in survey mode';
output.SWD_description          = 'RPW LFR L2s CWF science electric (potential difference) data in survey mode, time-tagged';
output.SWD_level                = 'L2S';
output.SWD_release_date         = '2016-06-02';
output.SWD_release_modification = INITIAL_RELEASE_MODIFICATION_STR;
sw_mode.outputs{end+1} = output;

C.sw_modes{end+1} = sw_mode;
%---------------------------------------------------------------------------------------------------
sw_mode = [];
%sw_mode.CLI_parameter = 'L2S_LFR-SURV-SWF-E';
sw_mode.CLI_parameter = 'LFR-SWF-E';
sw_mode.SWD_purpose = 'Generate snapshow waveform electric (potential difference) data from LFR';

sw_mode.inputs = {};

input = [];
input.CLI_parameter_name  = 'input_HK';
input.dataset_ID          = 'ROC-SGSE_HK_RPW-BIA';
input.dataset_version_str = '01';
sw_mode.inputs{end+1} = input;

input = [];
input.CLI_parameter_name  = 'input_SCI';
input.dataset_ID          = 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF';
input.dataset_version_str = '01';
sw_mode.inputs{end+1} = input;

sw_mode.outputs = {};

output = [];
output.JSON_output_file_identifier = 'output_SCI';
output.dataset_ID                  = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E';
output.dataset_version_str         = '01';
output.master_cdf_filename         = 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01.cdf';
output.SWD_name                 = 'LFR L2s SWF science electric data in survey mode';
output.SWD_description          = 'RPW LFR L2s SWF science electric (potential difference) data in survey mode, time-tagged';
output.SWD_level                = 'L2S';
output.SWD_release_date         = '2016-06-02';
output.SWD_release_modification = INITIAL_RELEASE_MODIFICATION_STR;
sw_mode.outputs{end+1} = output;

C.sw_modes{end+1} = sw_mode;
%---------------------------------------------------------------------------------------------------

constants = C;

end

